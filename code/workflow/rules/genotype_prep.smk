"""
snakemake -n -R geno_prep
snakemake --jobs 5 --use-singularity --use-conda -R geno_prep
snakemake --dag -R  geno_prep | dot -Tsvg > ../results/img/control/dag_geno_prep.svg

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
  --jn job_g.{name}.{jobid}.sh \
  -R geno_prep && mv job_g.* logs/
"""

GATK_REF = [ 'mirang', 'mirleo' ]
c_gatk = config[ 'sif_gatk' ]
c_qc = config[ 'sif_qc' ]
c_popgen = config[ 'sif_popgen' ]
c_sim = config[ 'sif_sim' ]

rule geno_prep:
    input:
      expand("../data/genomes/{species}_partitions.tsv", species = GATK_REF),
      expand("../results/checkpoints/index_ref_{species}.check", species = GATK_REF)

rule faidx_index:
    input:
      fa = "../data/genomes/{species}.fa.gz"
    output:
      fai = "../data/genomes/{species}.fa.gz.fai"
    resources:
      mem_mb=8192
    container: c_gatk
    shell:
      """
      samtools faidx {input.fa}
      """

rule check_genome:
    input:
      fai = expand("../data/genomes/{species}.fa.gz.fai", species = GATK_REF)
    output:
      partition = expand("../data/genomes/{species}_partitions.tsv", species = GATK_REF)
    conda: "r_tidy"
    log:
      "logs/r_partition_ref_genomes.log"
    shell:
      """
      Rscript R/partition_ref_genomes.R 2> {log} 1> {log}
      """

rule index_genomes:
  label "Q_def_index_c_gatk"
    input:
      fa = "../data/genomes/{species}.fa.gz",
      fai = "../data/genomes/{species}.fa.gz.fai"
    output: temp( touch("../results/checkpoints/index_ref_{species}.check") )
    resources:
      mem_mb=8192
    container: c_gatk
    shell:
      """
      bwa index -a bwtsw {input.fa}
      """
