"""
snakemake -n --rerun-triggers mtime -R all_gt_qc
snakemake --jobs 5 \
  --use-singularity --singularity-args "--bind $CDATA" \
  --use-conda --rerun-triggers mtime -R all_gt_qc
snakemake --rerun-triggers mtime --dag -R all_gt_qc | dot -Tsvg > ../results/img/control/dag_qc.svg

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
  --jn job_qc.{name}.{jobid}.sh \
  -R all_gt_qc && mv job_qc.* logs/
"""

rule all_gt_qc:
    message:
      """
      Check sequencing quality (`fastqc`), mapping coverage
      (`bamcov` and `bamtools`),
      inter-specific contamination (`fastq-screen`),
      allelic imbalance for intra-specific contamination.
      """
    input:
      by_sample_ln = expand( ["../results/checkpoints/fastqc/{sample_ln}.check",
                              "../results/qc/fastq_screen/{sample_ln}_fw_screen.txt"],
                              sample_ln = SAMPLES_LN ),
      by_sample = expand(["../results/qc/coverage/{sample_id}_on_{ref}_coverage.tsv.gz",
                          "../results/qc/coverage/masks/{sample_id}_on_{ref}_covmask.bed.gz",
                          "../results/qc/bamstats/{sample_id}_on_{ref}.bamstats",
                          "../results/het/{ref}_{sample_id}.csv"],
                         sample_id = SAMPLES, ref = GATK_REF[0]),
      by_ref = expand( "../results/img/qc/{ref}_{set}_het_stats.pdf",
                       ref = GATK_REF[0], set = ["all", "mirang", "mirleo"] )
      # GATK_REF[0] <- subset to mirang for now for disc-usage

# ---  sequencing qc --------------------------------
rule fastqc:
    input:
      fq_fw = lambda wc: "../data/raw_sequences/" + get_sample_info(wc, what = "file_fw"),
      fq_rv = lambda wc: "../data/raw_sequences/" + get_sample_info(wc, what = "file_rv")
    output: 
      check = touch( "../results/checkpoints/fastqc/{sample_ln}.check" )
    benchmark:
      'benchmark/qc/fastqc_{sample_ln}.tsv'
    resources:
      mem_mb=10240
    container: c_qc
    shell:
      """
      fastqc {input.fq_fw} {input.fq_rv} -o ../results/qc/fastqc/
      """

rule subsample_fastq:
    input:
      fq_fw = lambda wc: "../data/raw_sequences/" + get_sample_info(wc, what = "file_fw"),
      fq_rv = lambda wc: "../data/raw_sequences/" + get_sample_info(wc, what = "file_rv")
    output:
      fq_fw = temp( "../data/raw_sequences/subsets/{sample_ln}_fw.fq.gz" ),
      fq_rv = temp( "../data/raw_sequences/subsets/{sample_ln}_rv.fq.gz" )
    benchmark:
      'benchmark/qc/subsample_{sample_ln}.tsv'
    resources:
      mem_mb=8192
    container: c_qc
    shell:
      """
      seqtk sample -s 42 {input.fq_fw} 100000 | bgzip > {output.fq_fw}
      seqtk sample -s 42 {input.fq_rv} 100000 | bgzip > {output.fq_rv}
      """

rule fastq_screen:
    input:
      fq_fw = "../data/raw_sequences/subsets/{sample_ln}_fw.fq.gz",
      fq_rv = "../data/raw_sequences/subsets/{sample_ln}_rv.fq.gz"
    output:
      fw = "../results/qc/fastq_screen/{sample_ln}_fw_screen.txt",
      rv = "../results/qc/fastq_screen/{sample_ln}_rv_screen.txt"
    benchmark:
      'benchmark/qc/fq_screen_{sample_ln}.tsv'
    resources:
      mem_mb=8192,
      threads=7
    threads: 7
    container: c_qc
    shell:
      """
      fastq_screen --conf fastq_screen.conf --threads 7 --outdir ../results/qc/fastq_screen/ --aligner bowtie2 {input.fq_fw} {input.fq_rv}
      """

# ---  mapping qc -----------------------------------
rule merge_bam_by_sample:
    input:
      bams = lambda wc: "../results/mapped_bams/" + gather_sample_entries(wc, what = "sample_ln") + "_on_{ref}.dedup.bam"
    output:
      single_bam = '../results/mapped_bams/combined_per_sample/{sample_id}_on_{ref}.bam'
    benchmark:
      'benchmark/qc/merge_bam_{sample_id}_on_{ref}.tsv'
    resources:
      mem_mb=51200
    container: c_gatk
    shell:
      """
      BAMS=$( echo {input.bams} | sed "s/\[//g; s/\]//g; s/,//g; s/'//g" )
      samtools merge  -o {output.single_bam} $BAMS
      gatk --java-options "-Xmx45g" BuildBamIndex -I {output.single_bam} -O {output.single_bam}.bai
      """

rule bamcov:
    input:
      bam = "../results/mapped_bams/combined_per_sample/{sample_id}_on_{ref}.bam"
    output:
      cov = "../results/qc/coverage/{sample_id}_on_{ref}_coverage.tsv.gz"
    benchmark:
      'benchmark/qc/bamcov_{sample_id}_on_{ref}.tsv'
    params:
      unzipped = "../results/qc/coverage/{sample_id}_on_{ref}_coverage.tsv"
    resources:
      mem_mb=51200
    container: c_qc
    shell:
      """
      bamcov {input.bam} -o {params.unzipped}
      gzip {params.unzipped}
      """

rule bamcov_mask:
    input:
      bam = "../results/mapped_bams/combined_per_sample/{sample_id}_on_{ref}.bam"
    output:
      bed = "../results/qc/coverage/masks/{sample_id}_on_{ref}_covmask.bed.gz"
    benchmark:
      'benchmark/qc/covmask_{sample_id}_on_{ref}.tsv'
    resources:
      mem_mb=40960
    container: c_popgen
    shell:
      """
      samtools view -b {input.bam} | \
        genomeCoverageBed -ibam stdin -bg | \
        gzip > {output.bed}
      """

rule bamtools:
    input:
      bam = "../results/mapped_bams/combined_per_sample/{sample_id}_on_{ref}.bam"
    output:
      stats = "../results/qc/bamstats/{sample_id}_on_{ref}.bamstats"
    benchmark:
      'benchmark/qc/bamtools_{sample_id}_on_{ref}.tsv'
    resources:
      mem_mb=40960
    container: c_qc
    shell:
      """
      bamtools stats -in {input.bam} | \
        grep -v "*" > {output.stats}
      """

# ---  allelic imbalance ----------------------------
rule filter_mac1:
    input:
      vcf = "../results/genotyping/filtered/{ref}_bi-allelic.vcf.gz"
    output:
      vcf = "../results/genotyping/filtered/{ref}_mac1.vcf.gz"
    benchmark:
      'benchmark/qc/mac1_{ref}.tsv'
    resources:
      mem_mb=5120
    container: c_popgen
    shell:
      """
      vcftools \
        --gzvcf {input.vcf} \
        --max-mac 1 \
        --mac 1 \
        --recode \
        --stdout |  \
        bgzip > {output.vcf}
      tabix -p vcf {output.vcf}
      """

# rule all_pop:
#     input:
#       vcf = "../results/genotyping/filtered/{ref}_mac1.vcf.gz"
#     output:
#       inds = "../results/{ref}_all_inds.pop"
#     resources:
#       mem_mb=25600
#     container: c_popgen
#     shell:
#       """
#       # Creates individual file
#       zgrep "#CHROM" {input.vcf} | \
#         cut -f 10- | \
#         sed 's/\\t/\\n/g' > {output.inds}
#       """

rule export_het_ind:
    input:
      vcf = "../results/genotyping/filtered/{ref}_bi-allelic.vcf.gz",
      inds = "../results/inds_all.pop"
    output:
      hets = expand( "../results/het/{{ref}}_{sample_id}.csv", sample_id = SAMPLES )
    params:
      vcf_base = "{ref}_mac1",
      n_samples = len(SAMPLES),
      het_base = "../results/het/"
    benchmark:
      'benchmark/qc/het_{ref}.tsv'
    resources:
      mem_mb=25600
    container: c_popgen
    shell:
      """
      # Writes the AD fields from all heterozygous genotypes with minor allele count i into a file hets.i
      n={params.n_samples}
      c=1
      for i in $(cat {input.inds})
      do
        echo -ne "$i | $c of $n"
        vcftools \
          --gzvcf {input.vcf} \
          --mac 1 \
          --max-mac 1 \
          --indv $i \
           --recode \
           --stdout | \
            grep -v "^#" | \
             grep "0[/|]1" |\
              awk '{{split($9,a,":");for(i=1;i<=10;i++){{if(match(a[i],"AD")){{adidx=i}} }};for(i=10;i<=NF;i++){{if(match($i,"0[/|]1")==1){{split($i,b,":"); print b[adidx]}} }} }}' \
              > {params.het_base}/{wildcards.ref}_$i".csv"
        ((c=c+1))
      done
      """

rule bin_hets:
    input:
      inds = "../results/inds_{set}.pop",
      hets = expand( "../results/het/{{ref}}_{sample_id}.csv", sample_id = SAMPLES )
    output:
      freq = "../results/qc/allelic_imbalance/{ref,\w+}_{set,\w+}_het_ind_stats_freq2d.tsv",
      d = "../results/qc/allelic_imbalance/{ref,\w+}_{set,\w+}_het_ind_stats_d.tsv"
    params:
      het_base = "../results/het/{ref}_"
    benchmark:
      "benchmark/qc/allelic_imbalance_bin_{ref}_{set}.tsv"
    log:
      "logs/r_allelic_imbalance_bin_{ref}_{set}"
    conda: "r_tidy"
    resources:
      mem_mb=20480
    shell:
      """
      Rscript R/allelic_imbalance_het_bin.R \
        {params.het_base} \
        {input.inds} \
        {output.freq} \
        {output.d} 2> {log} 1> {log}
      """

rule plot_allelic_imbalance:
    input:
      freq = "../results/qc/allelic_imbalance/{ref}_{set}_het_ind_stats_freq2d.tsv",
      d = "../results/qc/allelic_imbalance/{ref}_{set}_het_ind_stats_d.tsv"
    output:
      plt = "../results/img/qc/{ref,\w+}_{set,\w+}_het_stats.pdf"
    benchmark:
      "benchmark/qc/allelic_imbalance_plot_{ref}_{set}.tsv"
    log:
      "logs/r_allelic_imbalance_plot_{ref}_{set}"
    conda: "r_tidy"
    resources:
      mem_mb=12288
    shell:
      """
      Rscript R/allelic_imbalance_plot.R \
        {input.freq} \
        {input.d} \
        {wildcards.ref} \
        {output.plt} 2> {log} 1> {log}
      """