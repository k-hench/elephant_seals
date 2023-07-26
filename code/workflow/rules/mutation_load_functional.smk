"""
snakemake --rerun-triggers mtime -n all_ml_snpeff
snakemake --dag  --rerun-triggers mtime -R all_ml_snpeff | dot -Tsvg > ../results/img/control/dag_snpeff_db.svg

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
  -R all_ml_snpeff
"""

GTF_FILE = "../data/genomes/annotation/mirang.gtf.gz"
GFF_FILE = "../data/genomes/annotation/mirang.gff3.gz"

rule all_ml_snpeff:
    input: 
      vcf = expand( "../results/genotyping/annotated/{vcf_pre}_ann.vcf.gz", vcf_pre = "mirang_filtered" )

rule download_gtf:
    output:
      gtf = GTF_FILE
    shell:
      """
      wget https://ftp.ncbi.nih.gov/genomes/refseq/vertebrate_mammalian/Mirounga_angustirostris/annotation_releases/current/GCF_021288785.2-RS_2023_03/GCF_021288785.2_ASM2128878v3_genomic.gtf.gz
      mv GCF_021288785.2_ASM2128878v3_genomic.gtf.gz {output.gtf}
      """

rule download_gff:
    output:
      gff = GFF_FILE
    shell:
      """
      wget https://ftp.ncbi.nih.gov/genomes/refseq/vertebrate_mammalian/Mirounga_angustirostris/annotation_releases/current/GCF_021288785.2-RS_2023_03/GCF_021288785.2_ASM2128878v3_genomic.gff.gz
      mv GCF_021288785.2_ASM2128878v3_genomic.gff.gz {output.gff}
      """

rule create_snpeff_config:
    output:
      conf = "../results/mutation_load/snp_eff/snpEff.config"
    shell:
      """
      echo "# Mirounga angustirostris genome, version GCF_021288785.2" > {output.conf}
      echo "mirang.genome : mirang" >> {output.conf}
      """

rule extract_cds:
    input:
      fa = "../data/genomes/mirang.fa",
      gff = "../data/genomes/annotation/mirang.gff3.gz",
    output:
      cds = "../results/mutation_load/snp_eff/data/mirang/cds.fa.gz"
    params:
      cds_prefix = "../results/mutation_load/snp_eff/data/mirang/"
    conda: "gff3toolkit"
    shell:
      """
      gff3_to_fasta \
        -g {input.gff} \
        -f {input.fa} \
        -st cds \
        -d complete \
        -o {params.cds_prefix}/mirang
      
      mv {params.cds_prefix}/mirang_cds.fa {params.cds_prefix}/cds.fa
      gzip {params.cds_prefix}/cds.fa
      """

rule extract_prot:
    input:
      fa = "../data/genomes/mirang.fa",
      gff = "../data/genomes/annotation/mirang.gff3.gz"
    output:
      pep = "../results/mutation_load/snp_eff/data/mirang/protein.fa.gz"
    params:
      pep_prefix = "../results/mutation_load/snp_eff/data/mirang"
    conda: "gff3toolkit"
    shell:
      """
      gff3_to_fasta \
        -g {input.gff} \
        -f {input.fa} \
        -st pep \
        -d complete \
        -o {params.pep_prefix}/mirang
      
      mv {params.pep_prefix}/mirang_pep.fa {params.pep_prefix}/protein.fa
      gzip {params.pep_prefix}/protein.fa
      """

rule create_snpeff_db:
    input:
      fa = "../data/genomes/mirang.fa",
      gtf = GTF_FILE,
      cds = "../results/mutation_load/snp_eff/data/mirang/cds.fa.gz",
      prot = "../results/mutation_load/snp_eff/data/mirang/protein.fa.gz",
      conf = "../results/mutation_load/snp_eff/snpEff.config"
    output:
      snp_fa = "../results/mutation_load/snp_eff/data/genomes/mirang.fa",
      snp_gff = "../results/mutation_load/snp_eff/data/mirang/genes.gtf.gz"
    params:
      snpeff_path = "../results/mutation_load/snp_eff"
    resources:
      mem_mb=25600
    container: c_ml
    shell:
      """
      mkdir -p {params.snpeff_path}/data/mirang {params.snpeff_path}/data/genomes
      cd {code_dir}/{params.snpeff_path}/data/mirang
      ln -s {code_dir}/{input.gtf} ./genes.gtf.gz
      cd {code_dir}/{params.snpeff_path}/data/genomes
      ln -s {code_dir}/{input.fa} ./mirang.fa
      cd {code_dir}/{params.snpeff_path}
      snpEff build -Xmx24G -c {code_dir}/{input.conf} -dataDir $(pwd)/data -gtf22 -v mirang
      """

rule run_snpeff:
    input:
      snpeff_gff = "../results/mutation_load/snp_eff/data/mirang/genes.gtf.gz",
      vcf = "../results/genotyping/filtered/{vcf_pre}.vcf.gz"
    output:
      snpef_vcf = "../results/genotyping/annotated/{vcf_pre}_ann.vcf",
      report = "../results/mutation_load/snp_eff/{vcf_pre}_stats.html"
    params:
      snpeff_path = "../results/mutation_load/snp_eff"
    resources:
      mem_mb=25600
    container: c_ml
    shell:
      """
      cd {code_dir}/{params.snpeff_path}
      ln -s {code_dir}/{input.vcf} ./
      snpEff ann -stats {wildcards.vcf_pre}_stats.html \
          -no-downstream \
          -no-intergenic \
          -no-intron \
          -no-upstream \
          -no-utr \
          -v \
          mirang {wildcards.vcf_pre}.vcf.gz > {code_dir}/{output.snpef_vcf}
      """

rule bgzip_vcf:
    input:
      vcf = "../results/genotyping/annotated/{vcf_pre}_ann.vcf"
    output:
      vcf = "../results/genotyping/annotated/{vcf_pre}_ann.vcf.gz",
      vcf_idx = "../results/genotyping/annotated/{vcf_pre}_ann.vcf.gz.tbi"
    conda: "popgen_basics"
    shell:
      """
      bgzip {input.vcf}
      tabix -p vcf {output.vcf}
      """