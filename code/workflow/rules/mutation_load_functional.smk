"""
snakemake --rerun-triggers mtime -n download_gtf

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
  --jn job_ml.{name}.{jobid}.sh \
  -R create_snpeff_db
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
      gtf = GTF_FILE,
      conf = "../results/mutation_load/snp_eff/snpEff.config"
    output:
      snp_fa = "../results/mutation_load/snp_eff/data/genomes/mirang.fa",
      snp_gff = "../results/mutation_load/snp_eff/data/mirang/genes.gtf.gz"
    params:
      snpeff_path = "../results/mutation_load/snp_eff"
    resources:
      mem_mb=15360
    container: c_ml
    shell:
      """
      mkdir -p {params.snpeff_path}/data/mirang {params.snpeff_path}/data/genomes
      cd {code_dir}/{params.snpeff_path}/data/mirang
      ln -s {code_dir}/{input.gtf} ./genes.gtf.gz
      cd {code_dir}/{params.snpeff_path}/data/genomes
      ln -s {code_dir}/{input.fa} ./mirang.fa
      cd {code_dir}/{params.snpeff_path}
      snpEff build  -c {code_dir}/{input.conf} -dataDir $(pwd) -gtf22 -v mirang
      """