"""
snakemake -n -R gt_all
snakemake --jobs 5 --use-singularity --singularity-args "--bind $CDATA" --use-conda -R gt_all
snakemake --dag -R  gt_all | dot -Tsvg > ../results/img/control/dag_genotyping.svg
snakemake --dag -R  gt_all --batch gt_all=1/30| dot -Tsvg > ../results/img/control/dag_genotyping_single.svg

snakemake --dag ../results/qc/snp_metrics/mirang_snp_metrics.tsv | dot -Tsvg > ../results/img/control/dag_mirang_metrics.svg

snakemake -n -R mapping_done
snakemake --jobs 5 --use-singularity --singularity-args "--bind $CDATA" --use-conda -R mapping_done
snakemake --dag -R  mapping_done | dot -Tsvg > ../results/img/control/dag_mapping.svg

snakemake --jobs 60 \
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
  -R mapping_done && mv job_gt.* logs/
"""

# read in the sequencing meta-data
seq_file_data = pd.read_table('../data/file_info.tsv')
seq_file_data['sample_ln'] = seq_file_data['sample_id'] + "_" +  [x[-1:] for x in seq_file_data['lane']]

SAMPLES_LN = seq_file_data['sample_ln'].values
SAMPLES = list(set(seq_file_data['sample_id'].values ))
GENOME_PARTITIONS = [ str(x + 1).zfill(2) for x in np.arange(20)]

# get a single entry from the file info table given a
# combination of sample_id + lane number
def get_sample_info(wildcards, what):
    return( seq_file_data[what][seq_file_data["sample_ln"] == wildcards.sample_ln].values[0] )

# get a set of all entries of a certain type (eg files)
# from the file info table given a sample_id
def gather_sample_enties(wildcards, what):
    return( seq_file_data[what][ seq_file_data["sample_id"] == wildcards.sample_id ].values )

rule gt_all:
    message:
      """
      producing the final filtered `vcf file`
      """
    input:
      vcf = expand("../results/genotyping/raw/{ref}_raw_snps.vcf.gz", ref = GATK_REF[0] )
      # GATK_REF[0] <- subset to mirang for now for disc-usage

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
    benchmark:
      'benchmark/genotyping/adapter_{sample_ln}.tsv'
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
    benchmark:
      'benchmark/genotyping/map_{sample_ln}_on_{ref}.tsv'
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
    benchmark:
      'benchmark/genotyping/merge_{sample_ln}_on_{ref}.tsv'
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
  benchmark:
    'benchmark/genotyping/dedup_{sample_ln}_on_{ref}.tsv'
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

rule collect_sample_and_haplotypecaller:
    input:
      bams = lambda wc: "../results/mapped_bams/" + gather_sample_enties(wc, what = "sample_ln") + "_on_{ref}.dedup.bam",
      ref = "../data/genomes/filtered/{ref}_filt.fa.gz"
    output:
      gvcf = "../results/genotypes/gvcf/{sample_id}_on_{ref}.g.vcf.gz"
    benchmark:
      'benchmark/genotyping/gvcf_{sample_id}_on_{ref}.tsv'
    resources:
      mem_mb=133120
    container: c_gatk
    shell:
      """
      INPUT="-I "$( echo {input.bams} | sed "s/\[//g; s/\]//g; s/,//g; s/'//g; s/ / -I /g" )
      gatk --java-options "-Xmx125g" \
        HaplotypeCaller  \
        -R {input.ref} \
        $INPUT \
        -O {output} \
        -ERC GVCF
      """

rule gather_gvcfs:
    input:
      gvcfs = expand( "../results/genotypes/gvcf/{sample_id}_on_{{ref}}.g.vcf.gz", sample_id = SAMPLES ),
      partitions = "../data/genomes/filtered/{ref}_filt_partitions.tsv"
    output:
      db =  directory( "../results/genotyping/{ref}_{part}_db" ),
      intervals = "../data/genomes/filtered/{ref}_{part}.intervals"
    benchmark:
      'benchmark/genotyping/gather_gvcf_{ref}_pt{part}.tsv'
    resources:
      mem_mb=133120
    container: c_gatk
    shell:
      """
      ALL_GVCF=$(echo " {input.gvcfs}" | sed "s/\[//g; s/\]//g; s/,//g; s/'//g; s/ / -V /g")

      grep "{wildcards.part}$" {input.partitions} | \
        cut -f 1 > {output.intervals}
    
      gatk --java-options "-Xmx125g" \
        GenomicsDBImport \
        $ALL_GVCF \
        --genomicsdb-workspace-path {output.db} \
        --tmp-dir tmp/ \
        -L {output.intervals}
      """

rule consolidate_gvcfs:
    input:
      db = "../results/genotyping/{ref}_{part}_db",
      ref = "../data/genomes/filtered/{ref}_filt.fa.gz"
    output:
      vcf = temp( "../results/genotyping/raw/{ref}_{part}_raw.vcf.gz" )
    benchmark:
      "benchmark/genotyping/raw_vcf_{ref}_pt{part}.tsv"
    resources:
      mem_mb=153600
    container: c_gatk
    shell:
      """
      gatk --java-options "-Xmx145g" \
        GenotypeGVCFs \
        -R {input.ref} \
        -V gendb://{input.db} \
        -O {output.vcf} \
        --tmp-dir tmp/
      """

rule gather_vcfs:
    input:
      vcfs = expand( "../results/genotyping/raw/{{ref}}_{part}_raw.vcf.gz", part = GENOME_PARTITIONS ),
      ref = "../data/genomes/filtered/{ref}_filt.fa.gz"
    output:
      vcf_list = temp( "tmp/{ref}_vcf_list.txt" ),
      vcf_all = "../results/genotyping/raw/{ref}_raw.vcf.gz",
      vcf_snps = "../results/genotyping/raw/{ref}_raw_snps.vcf.gz"
    benchmark:
      "benchmark/genotyping/raw_vcf_{ref}_collect.tsv"
    resources:
      mem_mb=133120
    container: c_gatk
    shell:
      """
      $(echo "{input.vcfs}" | sed "s/\[//g; s/\]//g; s/,//g; s/'//g; s/ /\n/g")
      echo *_raw.vcf.gz > {output.vcf_list}
    
      gatk GatherVcfsCloud \
        -I {output.vcf_list} \
        -O {output.vcf_all}
    
      tabix -p vcf {output.vcf_all}
    
      gatk --java-options "-Xmx95G" \
        SelectVariants \
        -R {input.ref} \
        -V {output.vcf_all} \
        --select-type-to-include SNP \
        -O {output.vcf_snps}
      """

rule snp_metrics:
    input:
      vcf = "../results/genotyping/raw/{ref}_raw_snps.vcf.gz"
    output:
      metrics = "../results/qc/snp_metrics/{ref}_snp_metrics.tsv"
    benchmark:
      "benchmark/genotyping/snp_metrics_{ref}.tsv"
    resources:
      mem_mb=30720
    container: c_gatk
    shell:
      """
      gatk --java-options "-Xmx25G" \
        VariantsToTable \
        --variant {input.vcf} \
        --output {output.metrics} \
        -F CHROM \
        -F POS \
        -F MQ \
        -F QD \
        -F FS \
        -F SOR \
        -F MQRankSum \
        -F ReadPosRankSum \
        --show-filtered
      """
