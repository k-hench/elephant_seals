"""
snakemake -n --rerun-triggers mtime -R geno_prep
snakemake --jobs 5 \
  --use-singularity --singularity-args "--bind $CDATA" \
  --use-conda --rerun-triggers mtime -R geno_prep
snakemake --dag -R  geno_prep | dot -Tsvg > ../results/img/control/dag_geno_prep.svg

snakemake --jobs 5 \
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
  --jn job_g.{name}.{jobid}.sh \
  -R geno_prep && mv job_g.* logs/
"""

GATK_REF = [ 'mirang', 'mirleo' ]
FASTQSCREEN_REF = deepcopy( GATK_REF )
FASTQSCREEN_REF.insert( 0, "galgal" )

rule geno_prep:
    input:
      expand("../data/genomes/{species}_partitions.tsv", species = GATK_REF),
      expand("../results/checkpoints/index_ref_{species}.check", species = GATK_REF),
      expand("../data/fq_screen_db/{species}", species = FASTQSCREEN_REF),
      expand("../data/fq_screen_db/{species}", species = FASTQSCREEN_REF),
      expand("../data/genomes/filtered/{species}_filt.dict", species = GATK_REF),
      '../results/checkpoints/prj_structure_ready.check'

rule faidx_index:
    input:
      fa = "../data/genomes/{species}.fa.gz"
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
      fai = expand("../data/genomes/{species}.fa.gz.fai", species = ['mirang', 'mirleo', 'filtered/mirang_filt', 'filtered/mirleo_filt'])
    output:
      partition = expand("../data/genomes/{species}_partitions.tsv", species = ['mirang', 'mirleo', 'filtered/mirang_filt', 'filtered/mirleo_filt'])
    container: c_conda
    conda: "r_tidy"
    log:
      "logs/r_partition_ref_genomes.log"
    shell:
      """
      Rscript R/partition_ref_genomes.R 2> {log} 1> {log}
      """

rule filter_genome:
    input: 
      fa = "../data/genomes/{species}.fa.gz",
    output:
      fa_filtered = "../data/genomes/filtered/{species,[a-z]+}_filt.fa.gz",
      fai_filtered = "../data/genomes/filtered/{species,[a-z]+}_filt.fa.gz.fai"
    resources:
      mem_mb=8192
    conda: 'bioawk'
    shell:
      """
      mkdir -p ../data/genomes/filtered/
      # bioawk filters by length and drops comments from fa headers
      bioawk -c fastx '{{ if(length($seq) > 1000) {{ print ">"$name; print $seq }} }}' {input.fa} | \
        fold | \
        bgzip > {output.fa_filtered}
      
      samtools faidx {output.fa_filtered}
      """

rule subset_genome_partitions:
    input: 
      partition = "../data/genomes/filtered/{species}_filt_partitions.tsv",
      fai_filtered = "../data/genomes/filtered/{species}_filt.fa.gz.fai"
    output:
      plt = "../results/img/qc/partition_sub_{species}.pdf",
      intervals = expand("../data/genomes/filtered/{{species}}_filt_partitions/part_{part}_sub_{sub}.intervals",
               part = GENOME_PARTITIONS, sub = (np.arange(10) + 1))
    params:
      n_partitions = 20,
      n_subs = 10
    container: c_conda
    conda: "r_tidy"
    log:
      "logs/r_partition_{species}_subset.log"
    shell:
      """
      Rscript R/partition_ref_genomes_subsets.R \
        {wildcards.species} \
        {params.n_partitions} \
        {params.n_subs} \
        {output.plt} 2> {log} 1> {log}
      """

rule bwa_index:
    input:
      fa = "../data/genomes/filtered/{species}_filt.fa.gz",
      fai = "../data/genomes/filtered/{species}_filt.fa.gz.fai"
    output: touch("../results/checkpoints/index_ref_{species}.check")
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

rule create_genome_dictionary:
  input:
    fa="../data/genomes/filtered/{species}_filt.fa.gz"
  output:
    dct="../data/genomes/filtered/{species}_filt.dict"
  resources:
    mem_mb=11264
  container: c_gatk
  shell:
    """
    gatk --java-options "-Xmx10G" \
      CreateSequenceDictionary \
      -R {input.fa} \
      -O {output.dct}
    """
