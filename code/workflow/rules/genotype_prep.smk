"""
snakemake -n -R geno_prep
snakemake --jobs 5 --use-singularity --singularity-args "--bind $CDATA" --use-conda -R geno_prep
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
FASTQSCREEN_REF = deepcopy( GATK_REF )
FASTQSCREEN_REF.insert( 0, "galgal" )

c_gatk = config[ 'sif_gatk' ]
c_qc = config[ 'sif_qc' ]
c_popgen = config[ 'sif_popgen' ]
c_sim = config[ 'sif_sim' ]

rule geno_prep:
    input:
      expand("../data/genomes/{species}_partitions.tsv", species = GATK_REF),
      expand("../results/checkpoints/index_ref_{species}.check", species = GATK_REF),
      expand("../data/fq_screen_db/{species}", species = FASTQSCREEN_REF),
      expand("../data/fq_screen_db/{species}", species = FASTQSCREEN_REF)

rule faidx_index:
    input:
      fa = "../data/genomes/{species,[a-z]+}.fa.gz"
    output:
      fai = "../data/genomes/{species,[a-z]+}.fa.gz.fai"
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


rule filter_genome:
    input: 
      fa = "../data/genomes/{species,[a-z]+}.fa",
      fai = "../data/genomes/{species,[a-z]+}.fa.gz.fai"
    output:
      fa_filtered = "../data/genomes/filtered/{species,[a-z]+}_filt.fa.gz",
      fai_filtered = "../data/genomes/filtered/{species,[a-z]+}_filt.fa.gz.fai",
      bed = '../results/genome/{species}_subset_500bp.bed'
    resources:
      mem_mb=8192
    container: c_popgen
    shell:
      """
      mkdir -p ../data/genomes/filtered/
      awk  -v OFS="\t" '$2 > 500 {{print \$1,0,\$2,\$1}}' {input.fai} > {output.bed}
      
      bedtools getfasta \
          -fi {input.fa} \
          -bed {output.bed} \
          -fullHeader \
          -nameOnly | \
          sed 's/=/ /g' | \
          bgzip > {output.fa_filtered}
      
      samtools faidx {output.fa_filtered}
      """

rule bwa_index:
    input:
      fa = "../data/genomes/filtered/{species}_filt.fa.gz",
      fai = "../data/genomes/filtered/{species}_filt.fa.gz.fai"
    output: temp( touch("../results/checkpoints/index_ref_{species}.check") )
    log:
      "logs/bwa_log/{species}.log"
    resources:
      mem_mb=8192
    container: c_gatk
    shell:
      """
      bwa index -a bwtsw {input.fa} &> {log}
      """

rule index_bowtie:
    input:
      fa = "../data/genomes/{species}.fa"
    output:
      bt_index = directory("../data/fq_screen_db/{species}/")
    params:
      bt_base = "../data/fq_screen_db/{species}/{species}"
    log:
      "logs/bt_log/{species}.log"
    resources:
      mem_mb=8192
    container: c_qc
    shell:
      """
      mkdir -p {output.bt_index}
      bowtie2-build {input} {params.bt_base} &> {log}
      """
