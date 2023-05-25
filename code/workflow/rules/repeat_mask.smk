'''
snakemake -n -R all_repeat

snakemake --jobs 10 \
  --latency-wait 30 \
  --use-singularity \
  --use-conda \
  -p \
  --configfile workflow/config.yml \
  --default-resources mem_mb=25600 \
  --cluster '
    qsub \
      -V -cwd \
      -P fair_share \
      -l idle=1 \
      -l si_flag=1 \
      -pe multislot {threads} \
      -l vf={resources.mem_mb}' \
      --jn job.{name}.{jobid}.sh \
      -R all_repeat && mv job.* logs/
'''

MASK_GENOMES=""

rule all_repeat:
    input: expand("../data/genomes/{genome}_hardmasked.fa.gz", genome = MASK_GENOMES )

rule build_db:
    input: '../data/genomes/{genome}.fa.gz'
    output: touch("../data/genomes/waypoints/db_{genome}.done")
    log:
      'logs/db/{genome}.log'
    conda: 'repeat_masking'
    container: None
    shell:
      """
      mkdir -p ../data/genomes/{wildcards.genome}_mod/ 
      BuildDatabase -name ../data/genomes/{wildcards.genome}_mod/{wildcards.genome}_db -engine ncbi ../data/genomes/{wildcards.genome}.fa.gz 2>{log}
      """

rule model_repeats:
    input:
      db_done = "../data/genomes/waypoints/db_{genome}.done"
    output:
      directory("../data/genomes/{genome}_mod/{genome}_model")
    params:
      wd = os.getcwd()
    log:
      'logs/mod/{genome}.log'
    conda: 'repeat_masking'
    container: None
    threads: 8
    shell:
      """
      cd ../data/genomes/{wildcards.genome}_mod/
      RepeatModeler -engine ncbi -pa {threads} -database {wildcards.genome}_db >& {params.wd}/{log}
      mv RM_* {wildcards.genome}_model
      """

rule mask_repeats:
    input: 
      model = "../data/genomes/{genome}_mod/{genome}_model",
      unmasked_genome = '../data/genomes/{genome}.fa.gz'
    output: '../data/genomes/{genome}_masked.fa.gz'
    log:
      'logs/mask/{genome}.log'
    conda: 'repeat_masking'
    container: None
    threads: 8
    shell:
      """
      mkdir -p ../data/genomes/tmp
      zcat {input.unmasked_genome} > ../data/genomes/tmp/{wildcards.genome}.tmp.fa
      RepeatMasker \
        -e rmblast \
        -pa {threads} -s \
        -lib {input.model}/consensi.fa.classified \
        -xsmall ../data/genomes/tmp/{wildcards.genome}.tmp.fa 2> {log}
      mv ../data/genomes/tmp/{wildcards.genome}.tmp.fa.masked data/genomes/{wildcards.genome}_masked.fa
      rm ../data/genomes/tmp/{wildcards.genome}.tmp.fa
      gzip ../data/genomes/{wildcards.genome}_masked.fa
      """

rule convert_to_hardmasked:
    input: '../data/genomes/{genome}_masked.fa.gz'
    output: '../data/genomes/{genome}_hardmasked.fa.gz'
    container: None
    shell:
      """
      zcat {input}| sed '/[>*]/!s/[atgcn]/N/g' | gzip > {output}
      """