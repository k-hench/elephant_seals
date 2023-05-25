"""
snakemake -c 1 --use-conda cactus_prep --configfile workflow/config.yml
"""
JOBSTORE_PATH='../results/cactus/jobStore.img'

rule cactus_prep:
    input: '../results/checkpoints/done_round_0.txt'

rule parse_cactus_config:
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

rule parse_cactus_steps:
    input: "../results/cactus/cactus_instructions.sh"
    output: "../results/checkpoints/done_round_0.txt"
    log: "logs/cactus/parse_instructions.log"
    conda: "r_base"
    shell:
      """
      Rscript --vanilla R/parse_cactus_jobs.R &> {log} && touch {output}
      """
