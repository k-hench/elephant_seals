"""
snakemake -c 1 --use-conda cactus_prep
snakemake --dag -R  cactus_prep | dot -Tsvg > ../results/img/control/dag_cactus_prep.svg
"""
JOBSTORE_PATH='../results/cactus/jobStore.img'

rule cactus_prep:
    input: '../results/checkpoints/done_round_0.check'

rule reformat_arcgaz:
    input: '../data/genomes/arcgaz_anc_h1.fa.gz'
    output: '../data/genomes/arcgaz.fa.gz'
    conda: "msa_align"
    shell:
      """
      zcat "^>" {input} | sed 's/=/./g; s/;/:/g' | bgzip > {output}
      """

rule parse_cactus_config:
    input: '../data/genomes/arcgaz.fa.gz'
    output: "../results/checkpoints/jobstore_setup.txt"
    log: "logs/cactus/parse_config.log"
    params:
      genomes = SPEC_ALL
    script: "../../py/sm_cactus_input.py"

rule jobstore_setup:
    input: "../results/checkpoints/jobstore_setup.txt"
    output: JOBSTORE_PATH
    params:
      [config['cactus_sif'], 'results/cactus/{name}.txt'.format(name = P_NAME)]
    log: "logs/cactus/jobstore_setup.log"
    script: "../../sh/sm_cactus_jobstore.sh"

rule stepwise_instructions:
    input: JOBSTORE_PATH
    output: "../results/cactus/cactus_instructions.sh"
    params: [ config['cactus_sif'], '../results/cactus/{name}.txt'.format(name = P_NAME), CACTUS_CORES ]
    log: "logs/cactus/instructions.log"
    script: "../../sh/sm_cactus_instructions.sh"

rule parse_cactus_step_1:
    input: "../results/cactus/cactus_instructions.sh"
    output: touch( "../results/checkpoints/parse-instructions_step_1.check" )
    log: "logs/cactus/parse_instructions_step1.log"
    conda: "r_base"
    shell:
      """
      Rscript --vanilla R/parse_cactus_jobs_step_1.R &> {log}
      """

rule parse_cactus_step_2:
    input: "../results/checkpoints/parse-instructions_step_1.check"
    output: touch( "../results/checkpoints/done_round_0.check" ),
            touch( "../results/checkpoints/done_round_0_j0_s0.check" )
    log: "logs/cactus/parse_instructions_setp2.log"
    conda: "r_base"
    shell:
      """
      Rscript --vanilla R/parse_cactus_jobs_step_2.R &> {log}
      """