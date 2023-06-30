"""
snakemake -n --rerun-triggers mtime -R gt_all
snakemake --jobs 5 \
  --use-singularity --singularity-args "--bind $CDATA" \
  --use-conda --rerun-triggers mtime -R gt_all
snakemake --dag  --rerun-triggers mtime -R gt_all | dot -Tsvg > ../results/img/control/dag_gt.svg

snakemake --dag ../results/img/control/snp_metrics_mirang.pdf | dot -Tsvg > ../results/img/control/dag_mirang_metrics.svg

snakemake -n --rerun-triggers mtime -R mapping_done
snakemake --jobs 5 \
  --use-singularity --singularity-args "--bind $CDATA" \
  --use-conda --rerun-triggers mtime -R mapping_done
snakemake --dag -R  mapping_done | dot -Tsvg > ../results/img/control/dag_mapping.svg

snakemake --jobs 60 \
  --latency-wait 30 \
  -p \
  --default-resources mem_mb=51200 threads=1 \
  --use-singularity \
  --singularity-args "--bind $CDATA" \
  --use-conda \
  --rerun-triggers mtime \
  --cluster '
    qsub \
      -V -cwd \
      -P fair_share \
      -l idle=1 \
      -l si_flag=1 \
      -pe multislot {threads} \
      -l vf={resources.mem_mb}' \
  --jn job_gt.{name}.{jobid}.sh \
  -R gt_all && mv job_gt.* logs/

files with ref subset to single species (GATK_REF[0]):
 - genotyping.smk
 - genotyping_qc.smk
"""

# read in the sequencing meta-data
seq_file_data = pd.read_table('../data/file_info.tsv')
seq_file_data['sample_ln'] = seq_file_data['sample_id'] + "_" +  [x[-1:] for x in seq_file_data['lane']]

SAMPLES_LN = seq_file_data['sample_ln'].values
SAMPLES = list(set(seq_file_data['sample_id'].values ))

# get a single entry from the file info table given a
# combination of sample_id + lane number
def get_sample_info(wildcards, what):
    return( seq_file_data[what][seq_file_data["sample_ln"] == wildcards.sample_ln].values[0] )

# get a set of all entries of a certain type (eg files)
# from the file info table given a sample_id
def gather_sample_entries(wildcards, what):
    return( seq_file_data[what][ seq_file_data["sample_id"] == wildcards.sample_id ].values )

rule gt_all:
    message:
      """
      producing the final filtered `vcf file`
      """
    input:
      vcf = expand("../results/genotyping/filtered/{ref}_bi-allelic.vcf.gz", ref = GATK_REF[0] )
      # GATK_REF[0] <- subset to mirang for now for disc-usage

rule gt_invariant:
    message:
      """
      producing the final filtered `vcf file`
      """
    input:
      vcf = expand( "../results/genotyping/filtered/partitions/{ref}_all_bp_{part}_filtered.vcf.gz", part = GENOME_PARTITIONS, ref = GATK_REF[0] )
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
      bams = lambda wc: "../results/mapped_bams/" + gather_sample_entries(wc, what = "sample_ln") + "_on_{ref}.dedup.bam",
      ref = "../data/genomes/filtered/{ref}_filt.fa.gz"
    output:
      gvcf = "../results/genotyping/gvcf/{sample_id}_on_{ref,[^0-9]+}.g.vcf.gz"
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
      gvcfs = expand( "../results/genotyping/gvcf/{sample_id}_on_{{ref}}.g.vcf.gz", sample_id = SAMPLES ),
      partitions = "../data/genomes/filtered/{ref}_filt_partitions.tsv"
    output:
      db =  directory( "../results/genotyping/{ref,[^0-9]+}_{part,[0-9]*}_db" ),
      intervals = "../data/genomes/filtered/{ref,[^0-9]+}_{part,[0-9]*}.intervals"
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
      vcf = temp( "../results/genotyping/raw/{ref,[^0-9]+}_{part}_raw.vcf.gz" )
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
      vcf_list = temp( "tmp/{ref}_vcf.list" ),
      vcf_all = protected( "../results/genotyping/raw/{ref}_raw.vcf.gz" ),
      vcf_snps = "../results/genotyping/raw/{ref}_raw_snps.vcf.gz"
    benchmark:
      "benchmark/genotyping/raw_vcf_{ref}_collect.tsv"
    resources:
      mem_mb=133120
    container: c_gatk
    shell:
      """
      echo "{input.vcfs}" | sed "s/ /\\n/g" > {output.vcf_list}
    
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

rule metrics_density:
    input:
      metrics = "../results/qc/snp_metrics/{ref}_snp_metrics.tsv"
    output: touch("../results/checkpoints/gatk_snp_metrics_density_{ref}.check")
    benchmark:
      "benchmark/genotyping/snp_metrics_density_{ref}.tsv"
    log:
      "logs/r_gatk_metrics_density_{ref}.log"
    conda: "r_tidy"
    resources:
      mem_mb=20480
    shell:
      """
      Rscript R/gatk_metrics_density.R {wildcards.ref} {input.metrics} 2> {log} 1> {log}
      """

rule metrics_plot:
    input:
      check = "../results/checkpoints/gatk_snp_metrics_density_{ref}.check"
    output:
      plot = "../results/img/control/snp_metrics_{ref}.pdf"
    benchmark:
      "benchmark/genotyping/snp_metrics_density_plot_{ref}.tsv"
    log:
      "logs/r_gatk_metrics_plot_{ref}.log"
    conda: "r_tidy"
    resources:
      mem_mb=5120 
    shell:
      """
      Rscript R/gatk_metrics_plot.R {wildcards.ref} 2> {log} 1> {log}
      """

# data farme with the filter thresholds
# (currently dummy values 1-16, to be replaced once
# the snp metrics have completed)
filter_params = pd.DataFrame({'mirang': [7.5, 17.5, 55.0, 3.0, -0.5, 0.5, -2.25, 2.25],
                              'mirleo': [7.5, 17.5, 55.0, 3.0, -0.5, 0.5, -2.25, 2.25]},
                              index = ["qd", "fs", "mq", "sor",
                                       "mq_r_lower", "mq_r_upper",
                                       "rpos_lower", "rpos_upper"])

def get_filter_params(wildcards):
    return(filter_params[wildcards.ref])

rule gatk_filter_snps:
    input:
      vcf = "../results/genotyping/raw/{ref}_raw_snps.vcf.gz",
      ref = "../data/genomes/filtered/{ref}_filt.fa.gz",
      metrics_plot = "../results/img/control/snp_metrics_{ref}.pdf"
    output:
      vcf_flagged = temp( "../results/genotyping/raw/{ref}_flagged.vcf.gz" ),
      vcf_filtered = "../results/genotyping/filtered/{ref}_filtered.vcf.gz"
    params:
      vals = lambda wc: get_filter_params(wc)
    benchmark:
      "benchmark/genotyping/snp_filtering_{ref}.tsv"
    resources:
      mem_mb=51200
    container: c_gatk
    shell:
      """
      gatk --java-options "-Xmx35G" \
        VariantFiltration \
        -R {input.ref} \
        -V {input.vcf} \
        -O {output.vcf_flagged} \
        --filter-expression "QD < {params.vals[qd]}" \
        --filter-name "filter_QD" \
        --filter-expression "FS > {params.vals[fs]}" \
        --filter-name "filter_FS" \
        --filter-expression "MQ < {params.vals[mq]}" \
        --filter-name "filter_MQ" \
        --filter-expression "SOR > {params.vals[sor]}" \
        --filter-name "filter_SOR" \
        --filter-expression "MQRankSum < {params.vals[mq_r_lower]} || MQRankSum > {params.vals[mq_r_upper]}" \
        --filter-name "filter_MQRankSum" \
        --filter-expression "ReadPosRankSum < {params.vals[rpos_lower]} || ReadPosRankSum > {params.vals[rpos_upper]}" \
        --filter-name "filter_ReadPosRankSum"
    
      gatk --java-options "-Xmx35G" \
        SelectVariants \
        -R {input.ref} \
        -V {output.vcf_flagged} \
        -O {output.vcf_filtered} \
        --exclude-filtered
      """

rule vcftools_snp_filter:
    input:
      vcf = "../results/genotyping/filtered/{ref}_filtered.vcf.gz"
    output:
      vcf = "../results/genotyping/filtered/{ref}_bi-allelic.vcf.gz"
    benchmark:
      "benchmark/genotyping/snp_filtering_vcftools_{ref}.tsv"
    resources:
      mem_mb=40960
    container: c_popgen
    shell:
      """
      vcftools \
        --gzvcf {input.vcf} \
        --max-missing-count 5 \
        --max-alleles 2 \
        --stdout  \
        --recode | \
        bgzip > {output.vcf}
    
      tabix -p vcf {output.vcf}
      """

# ---  invariant site calling (all_bp) --------------
rule consolidate_gather_all_bp:
    input:
      db = "../results/genotyping/{ref}_{part}_db",
      ref = "../data/genomes/filtered/{ref}_filt.fa.gz",
      intervals = "../data/genomes/filtered/{ref}_filt_partitions/part_{part}_sub_{sub}.intervals"
    output:
      vcf_raw = temp( "../results/genotyping/all_bp/raw/{ref,[^0-9]+}_all_bp_{part,[0-9]*}_sub_{sub,[0-9]*}_raw.vcf.gz" )
    benchmark:
      'benchmark/genotyping/gather_gvcf_all_bp_{ref}_pt{part}_s{sub}.tsv'
    log:
      "logs/genotyping/all_bp/gather_{ref}_pt{part}_s{sub}.log"
    resources:
      mem_mb=174080
    container: c_gatk
    shell:
      """
      gatk --java-options "-Xmx135g" \
        GenotypeGVCFs \
        -R {input.ref} \
        -V gendb://{input.db} \
        -O {output.vcf_raw} \
        --tmp-dir tmp/ \
        -L {input.intervals} \
        --include-non-variant-sites true 2> {log} 1> {log}
      """

rule select_var_all_bp:
    input:
      ref = "../data/genomes/filtered/{ref}_filt.fa.gz",
      intervals = "../data/genomes/filtered/{ref}_filt_partitions/part_{part}_sub_{sub}.intervals",
      vcf_raw = "../results/genotyping/all_bp/raw/{ref}_all_bp_{part}_sub_{sub}_raw.vcf.gz"
    output:
      vcf_snps = "../results/genotyping/all_bp/raw/{ref,[^0-9]+}_all_bp_{part,[0-9]*}_sub_{sub,[0-9]*}_raw_snps.vcf.gz"
    benchmark:
      'benchmark/genotyping/select_var_all_bp_{ref}_pt{part}_sub_{sub}.tsv'
    resources:
      mem_mb=174080
    container: c_gatk
    shell:
      """
      gatk --java-options "-Xmx135g" \
        SelectVariants \
        -L {input.intervals} \
        -R {input.ref} \
        -V {input.vcf_raw} \
        --select-type-to-exclude INDEL \
        -O {output.vcf_snps} 
      """

rule gatk_filter_snps_all_bp:
    input:
      vcf = "../results/genotyping/all_bp/raw/{ref}_all_bp_{part}_sub_{sub}_raw_snps.vcf.gz",
      ref = "../data/genomes/filtered/{ref}_filt.fa.gz",
      metrics_plot = "../results/img/control/snp_metrics_{ref}.pdf"
    output:
      vcf_flagged = temp( "../results/genotyping/raw/{ref,[^0-9]+}_all_bp_{part,[0-9]*}_sub_{sub,[0-9]*}_flagged.vcf.gz" ),
      vcf_filtered = temp( "../results/genotyping/filtered/partitions/{ref,[^0-9]+}_all_bp_{part,[0-9]*}_sub_{sub,[0-9]*}_filtered.vcf.gz")
    params:
      vals = lambda wc: get_filter_params(wc)
    benchmark:
      'benchmark/genotyping/snp_filtering_all_bp_{ref}_pt{part}_s{sub}.tsv'
    resources:
      mem_mb=61440
    container: c_gatk
    shell:
      """
      gatk --java-options "-Xmx55G" \
        VariantFiltration \
        -R {input.ref} \
        -V {input.vcf} \
        -O {output.vcf_flagged} \
        --filter-expression "QD < {params.vals[qd]}" \
        --filter-name "filter_QD" \
        --filter-expression "FS > {params.vals[fs]}" \
        --filter-name "filter_FS" \
        --filter-expression "MQ < {params.vals[mq]}" \
        --filter-name "filter_MQ" \
        --filter-expression "SOR > {params.vals[sor]}" \
        --filter-name "filter_SOR" \
        --filter-expression "MQRankSum < {params.vals[mq_r_lower]} || MQRankSum > {params.vals[mq_r_upper]}" \
        --filter-name "filter_MQRankSum" \
        --filter-expression "ReadPosRankSum < {params.vals[rpos_lower]} || ReadPosRankSum > {params.vals[rpos_upper]}" \
        --filter-name "filter_ReadPosRankSum"
    
      gatk --java-options "-Xmx55G" \
        SelectVariants \
        -R {input.ref} \
        -V {output.vcf_flagged} \
        -O {output.vcf_filtered} \
        --exclude-filtered
      """

rule concat_parts_snps_all_bp:
    input: 
      vcfs = expand( "../results/genotyping/filtered/partitions/{{ref}}_all_bp_{{part}}_sub_{sub}_filtered.vcf.gz", sub = (np.arange(10) + 1) )
    output:
      vcf_list = temp( "../results/genotyping/filtered/partitions/{ref}_all_bp_{part}_filtered_vcf.list" ),
      vcf = "../results/genotyping/filtered/partitions/{ref}_all_bp_{part}_filtered.vcf.gz"
    shell:
      """
      echo "{input.vcfs}" | sed "s/ /\\n/g" > {output.vcf_list}
    
      gatk --java-options "-Xmx55G" \
        MergeVcfs \
        -I {output.vcf_list} \
        -O {output.vcf}
      """