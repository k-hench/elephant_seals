"""
snakemake -n -R geno_prep
snakemake -c 1 --use-singularity  -R geno_prep
snakemake --dag -R  geno_prep | dot -Tsvg > ../results/img/control/dag_geno_prep.svg
"""

GATK_REF = [ 'mirang', 'mirleo' ]
c_gatk = config[ 'sif_gatk' ]
c_qc = config[ 'sif_qc' ]
c_popgen = config[ 'sif_popgen' ]
c_sim = config[ 'sif_sim' ]

rule geno_prep:
    input: expand("data/genomes/{species}_partitions.tsv", species = GATK_REF)

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
      partition = expand("data/genomes/{species}_partitions.tsv", species = GATK_REF)
    conda: "r_tidy"
    log:
      "logs/r_partition_ref_genomes.log"
    shell:
      """
      Rscript R/partition_ref_genomes.R 2> {log} 1> {log}
      """
