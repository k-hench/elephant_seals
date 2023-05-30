"""
snakemake -n -R gt_all
snakemake --jobs 5 --use-singularity --singularity-args "--bind $CDATA" --use-conda -R gt_all
snakemake --dag -R  gt_all | dot -Tsvg > ../results/img/control/dag_genotyping.svg

snakemake -n -R mapping_done
snakemake --jobs 5 --use-singularity --singularity-args "--bind $CDATA" --use-conda -R mapping_done
snakemake --dag -R  mapping_done | dot -Tsvg > ../results/img/control/dag_mapping.svg

snakemake --jobs 5 \
  --latency-wait 30 \
  -p \
  --default-resources mem_mb=51200 threads=1 \
  --use-singularity \
  --singularity-args "--bind $CDATA" \
  --use-conda \
  --cluster '
    qsub \
      -V -cwd \
      -P fair_share \
      -l idle=1 \
      -l si_flag=1 \
      -pe multislot {threads} \
      -l vf={resources.mem_mb}' \
  --jn job_gt.{name}.{jobid}.sh \
  -R mapping_done && mv job_g.* logs/
"""

# read in the sequencing meta-data
seq_file_data = pd.read_table('../data/file_info.tsv')
seq_file_data = seq_file_data.set_index(seq_file_data['sample_id'] + "_" +  [x[-1:] for x in seq_file_data['lane']], drop = False)

SAMPLES_LN = seq_file_data.index.values

def get_sample_info(wildcards, what):
    return( seq_file_data.loc[wildcards.sample_ln, what ] )

rule gt_all:
    message:
      """
      producing the final filtered `vcf file`
      """

rule mapping_done:
    message:
      """
      All data has been mapped and de-duplicated -> ready for haplotypecaller.
      """
    input:
      deduped = expand( "../results/mapped_bams/{sample_ln}_on_{ref}.dedup.bam", sample_ln = SAMPLES_LN, ref = GATK_REF[0] )
      # GATK_REF[0] <- subset to mirang for now for disc-usage

# ---  actual genotyping ----------------------------

rule ubam_adapters:
    input:
      fq_fw = lambda wc: "../data/raw_sequences/" + get_sample_info(wc, what = "file_fw"),
      fq_rv = lambda wc: "../data/raw_sequences/" + get_sample_info(wc, what = "file_rv")
    output:
      ubam = temp( 'tmp/ubams/{sample_ln}.ubam.bam' ),
      adapter_bam = temp( 'tmp/ubams/{sample_ln}.adapter.bam' ), 
      adapter_fq = temp( 'tmp/fq/{sample_ln}.adapter.fq' )
    params:
      sample_id = lambda wc: get_sample_info(wc, what = "sample_id"),
      flowcell = lambda wc: get_sample_info(wc, what = "flowcell_id"),
      lane = lambda wc: get_sample_info(wc, what = "lane_header"),
      company = lambda wc: get_sample_info(wc, what = "company")
    resources:
      mem_mb=15360
    container: c_gatk
    shell:
      """
      # fastq to ubam
      gatk --java-options "-Xmx14G" \
        FastqToSam \
        --SAMPLE_NAME {params.sample_id} \
        -F1 {input.fq_fw} \
        -F2 {input.fq_rv} \
        -O {output.ubam} \
        -RG {wildcards.sample_ln} \
        -LB {params.sample_id}_lib1 \
        -PU {params.flowcell}_{params.lane} \
        -PL Illumina \
        -CN {params.company} \
        --TMP_DIR tmp/;
      
      # mark adapters
      gatk --java-options "-Xmx14G" \
        MarkIlluminaAdapters \
        -I {output.ubam} \
        -O {output.adapter_bam} \
        --M ../results/qc/adapter/{wildcards.sample_ln}.adapter.metrics.tsv \
        --TMP_DIR tmp/
      
      gatk --java-options "-Xmx14G" \
        SamToFastq \
        -I {output.adapter_bam}  \
        --FASTQ {output.adapter_fq}\
        --CLIPPING_ATTRIBUTE XT \
        --CLIPPING_ACTION 2 \
        --INTERLEAVE true \
        --NON_PF true \
        --TMP_DIR tmp/
      """

rule bwa_mapping:
    input:
      fq = 'tmp/fq/{sample_ln}.adapter.fq',
      ref = "../data/genomes/filtered/{ref}_filt.fa.gz"
    output:
      bam = temp( 'tmp/bam/{sample_ln}_on_{ref}.bam' )
    resources:
      mem_mb=15360
    threads: 5
    container: c_gatk
    shell:
      """
      source sh/params_switch.sh
      bwa mem \
        -M -t ${{CORES}} \
        -p {input.ref} \
        {input.fq} | \
        samtools view -bS - > {output.bam}
      """

rule merge_bam:
    input: 
      bam = 'tmp/bam/{sample_ln}_on_{ref}.bam',
      ubam = 'tmp/ubams/{sample_ln}.ubam.bam',
      ref = "../data/genomes/filtered/{ref}_filt.fa.gz"
    output:
      bam = temp( "tmp/merged/{sample_ln}_on_{ref}.bam" )
    resources:
      mem_mb=81920
    container: c_gatk
    shell:
      """
      gatk --java-options "-Xmx75G" \
        MergeBamAlignment \
        --ALIGNED_BAM {input.bam} \
        --UNMAPPED_BAM {input.ubam} \
        --OUTPUT {output.bam} \
        -R {input.ref} \
        --CREATE_INDEX true \
        --ADD_MATE_CIGAR true \
        --CLIP_ADAPTERS false \
        --CLIP_OVERLAPPING_READS true \
        --INCLUDE_SECONDARY_ALIGNMENTS true \
        --MAX_INSERTIONS_OR_DELETIONS -1 \
        --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
        --ATTRIBUTES_TO_RETAIN XS \
        --TMP_DIR tmp/
      """

rule mark_duplicates:
      input:
        bam = "tmp/merged/{sample_ln}_on_{ref}.bam",
        ref = "../data/genomes/filtered/{ref}_filt.fa.gz"
      output:
        marked_bam = temp( "tmp/dup_marked/{sample_ln}_on_{ref}.marked.bam" ),
        sorted_bam = temp( "tmp/sorted/{sample_ln}_on_{ref}.sorted.bam" ),
        final_bam = "../results/mapped_bams/{sample_ln}_on_{ref}.dedup.bam",
        metrics =  '../results/qc/dedup/{sample_ln}_on_{ref}_dedup_metrics.tsv'
      resources:
        mem_mb=122880
      container: c_gatk
      shell:
        """
        gatk --java-options "-Xmx110G" \
          MarkDuplicates \
          -I {input.bam} \
          -O {output.marked_bam} \
          -M {output.metrics} \
          -MAX_FILE_HANDLES 1000 \
          -TMP_DIR tmp/
        
        gatk --java-options "-Xmx110G" \
          SortSam \
          -I {output.marked_bam} \
          -O {output.sorted_bam} \
          --SORT_ORDER "coordinate" \
          --CREATE_INDEX false \
          --CREATE_MD5_FILE false \
          --TMP_DIR tmp/
        
        gatk --java-options "-Xmx110G" \
          SetNmMdAndUqTags \
          --INPUT {output.sorted_bam} \
          --OUTPUT {output.final_bam} \
          --CREATE_INDEX true \
          --TMP_DIR tmp/ \
          --REFERENCE_SEQUENCE {input.ref}
        """