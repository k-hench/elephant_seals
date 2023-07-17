"""
snakemake --rerun-triggers mtime -n download_gtf
"""

GTF_FILE = "../data/genomes/annotation/mirang.gtf.gz"

rule all_ml_snpeff:
    input: ""

rule download_gtf:
    output:
      gtf = GTF_FILE
    shell:
      """
      wget https://ftp.ncbi.nih.gov/genomes/refseq/vertebrate_mammalian/Mirounga_angustirostris/annotation_releases/current/GCF_021288785.2-RS_2023_03/GCF_021288785.2_ASM2128878v3_genomic.gtf.gz
      mv GCF_021288785.2_ASM2128878v3_genomic.gtf.gz {output.gtf}
      """

rule create_snpeff_config:
    output:
      conf = "../results/mutation_load/snp_eff/snpEff.config"
    shell:
      """
      echo "# Mirounga angustirostris genome, version GCF_021288785.2" > {output.conf}
      echo "mirang.genome : mirang" >> {output.conf}
      """

rule create_snpeff_db:
    input:
      fa = "../data/genomes/mirang.fa",
      gtf = GTF_FILE
      conf = "../results/mutation_load/snp_eff/snpEff.config"
    params:
      snpeff_path = "../results/mutation_load/snp_eff"
    resources:
      mem_mb=15360
    container: c_ml
    shell:
      """
      mkdir -p {params.snpeff_path}/data/mirang {params.snpeff_path}/data/genomes
      cd {params.snpeff_path}/data/mirang
      ln -s {input.gtf} ./genes.gtf.gz
      cd {params.snpeff_path}/data/genomes
      ln -s {input.fa} ./mirang.fa
      cd {params.snpeff_path}
      snpEff build  -c {input.conf} -dataDir $(pwd) -gtf22 -v mirang
      """