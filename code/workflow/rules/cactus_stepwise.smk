"""
snakemake --rerun-triggers mtime -n -R cactus_stepwise
snakemake --dag -R  cactus_stepwise | dot -Tsvg > ../results/img/control/dag_cactus_step.svg
snakemake --jobs 3 -R  cactus_stepwise

snakemake --jobs 50 \
  --latency-wait 30 \
  -p \
  --default-resources mem_mb=51200 threads=1 \
  --use-singularity \
  --singularity-args "--bind $CDATA" \
  --use-conda \
  --rerun-triggers mtime \
  --cluster '
    sbatch \
      --export ALL \
      --ntasks 1 \
      --cpus-per-task {threads} \
      --error logs/{name}.{jobid}.err \
      --output logs/{name}.{jobid}.out \
      --mem {resources.mem_mb} \
      --job-name {name}.{jobid}' \
      --jn job_c.{name}.{jobid}.sh \
      -R cactus_stepwise
"""
localrules: cactus_stepwise, round_completed, cactus_export_hal


job_file = "../results/cactus/job_inventory.tsv"
rounds = pd.read_table(job_file)['round']
n_rounds = rounds.max()
n_jobs = pd.read_table(job_file)['n_jobs']
s_bind_paths = config[ 'singularity_bind_paths' ]

job_tbl = pd.read_table("../results/cactus/job_list.tsv")
final_step_tbl = pd.read_table("../results/cactus/job_final_steps.tsv")

rule cactus_stepwise:
    input: '../results/cactus/{name}_check.check'.format(name = P_NAME)

def collect_steps(rnd, job):
    step = (final_step_tbl['step_idx'][(final_step_tbl["round_idx"] == rnd) & (final_step_tbl["job_idx"] == job)]).values[0]
    prev_steps = ((job_tbl["step_idx"][(job_tbl["round_idx"] == rnd) & (job_tbl["job_idx"] == job)]).reset_index(drop = True))[:step].values
    step_checks = [ '../results/checkpoints/done_round_' + str(rnd) + "_j" + str(job) + "_s" + str(i) + ".check" for i in prev_steps ]
    return( step_checks )

def collect_jobs(wildcards):
  rnd = int(wildcards.nr)
  if rnd == 0:
    return( "../results/checkpoints/done_round_0.check" )
  else:
    n_j = n_jobs[rnd - 1]
    j_list = (np.arange(0, n_j) + 1)
    steps = np.array([])
    for j in j_list:
      steps = np.append(steps, collect_steps(rnd, j))
    # j_checks = [ 'results/checkpoints/done_round_' + str(rnd) + "_j" + str(i) + ".check" for i in j_list ]
    return(steps)

def previous_round(wildcards):
  rnd = int(wildcards.nr)
  return("../results/checkpoints/done_round_" + str(rnd-1) + ".check")

# def steps_ancestors(wildcards):
#   step_list = (np.arange(0, n_j) + 1)

# def previous_step(wildcards):
#   rnd = int(wildcards.nr)
#   job = int(wildcards.job)
#   step = int(wildcards.step)
#   job_tbl[job_tbl['round_idx'] == rnd & job_tbl['job_idx'] == job]

rule round_completed:
    input: lambda wc: collect_jobs(wc)
    output: "../results/checkpoints/done_round_{nr}.check"
    shell:
      '''
      touch {output}
      '''

# parse_job(obj, "sh")
# obj = type('Test', (object,), {})
# obj.nr = 2
# obj.job = 1
# obj.step = 2
def parse_job(wildcards, output):
  # >> parse the job path <<
  # It is not genreally clear how the round and 
  # job fromat looks like: they are numbered with 
  # padding zeros depending on the overall number
  # of rounds/jobs.
  # The first round might therefore be for example
  # either "round_1" or round "round_001".
  rnd = int(wildcards.nr)
  job = int(wildcards.job)
  step = int(wildcards.step)
  cur_rnd = ((job_tbl["round"][job_tbl["round_idx"] == rnd]).reset_index(drop = True))[0]
  cur_job = ((job_tbl["job"][(job_tbl["round_idx"] == rnd) & (job_tbl["job_idx"] == job)]).reset_index(drop = True))[step - 1][:-3]
  cur_step = ((job_tbl["step_idx"][(job_tbl["round_idx"] == rnd) & (job_tbl["job_idx"] == job)]).reset_index(drop = True))[step - 1]
  cur_step_sh = ((job_tbl["step"][(job_tbl["round_idx"] == rnd) & (job_tbl["job_idx"] == job)]).reset_index(drop = True))[step - 1]
  prev_steps = ((job_tbl["step_idx"][(job_tbl["round_idx"] == rnd) & (job_tbl["job_idx"] == job)]).reset_index(drop = True))[:step-1].values
  j_checks = [ '../results/checkpoints/done_round_' + str(rnd) + "_j" + str(job) + "_s" + str(i) + ".check" for i in prev_steps ]
  if output == "sh":
    return( "sh/cactus/" + cur_rnd + "/" + cur_job + "/" + cur_step_sh + ".sh" )
  elif output == "prev":
    return( j_checks )
  elif output == "check":
    return( '../results/checkpoints/done_round_' + str(rnd) + "_j" + str(job) + "_s" + str(cur_step) + ".check" )
    # parse_job(pd.Series([2, 1], index = ["nr", "job"]))

rule single_job:
    input:
      previous_round = lambda wc: previous_round(wc),
      previous_steps = lambda wc: parse_job(wc, "prev"),
      job_script = lambda wc: parse_job(wc, "sh")
    output: "../results/checkpoints/done_round_{nr}_j{job}_s{step}.check"
    params:
      sif = config['cactus_sif'],
      seqfile = '../results/cactus/{name}.txt'.format(name = P_NAME),
      jobstore = JOBSTORE_PATH
    log: "logs/cactus/jobs/round_{nr}_j{job}_s{step}.log"
    threads: 24
    shell:
      '''
      readonly CACTUS_IMAGE={params.sif} 
      readonly SEQFILE={params.seqfile}
      readonly SEQNAME=${{SEQFILE##*/}}
      readonly RUN_ID=${{SEQNAME%.txt}}
      readonly CACTUS_SCRATCH=../results/cactus/scratch/${{RUN_ID}}

      echo "file: " ${{SEQFILE}} &>> {log}
      echo "==================" &>> {log}
      echo "img: "${{CACTUS_IMAGE}} &>> {log}
      echo "==================" &>> {log}
      echo "round {wildcards.nr}; job {wildcards.job}" &>> {log}
      echo "==================" &>> {log}
      echo "bind {s_bind_paths}" &>> {log}
      echo "==================" &>> {log}

      apptainer exec --cleanenv \
        --fakeroot --overlay ${{CACTUS_SCRATCH}} \
        --bind ${{CACTUS_SCRATCH}}/tmp:/tmp,$(pwd),{s_bind_paths} \
        --env PYTHONNOUSERSITE=1 \
        {params.sif} \
        bash {input.job_script} &>> {log}

      touch {output}
      '''

rule cactus_export_hal:
    input:
      final_step_check = "../results/checkpoints/done_round_{nr}.check".format(nr = n_rounds)
    output:
      hal = "../results/cactus/{name}.hal".format(name = P_NAME)
    params:
      hal = "../results/cactus/scratch/{name}/tmp/steps-output/{name}.hal".format(name = P_NAME)
    shell:
      """
      mv {params.hal} {output.hal}
      """

rule cactus_check:
    input: '../results/cactus/{name}.hal'.format(name = P_NAME)
    output: '../results/cactus/{name}_check.check'.format(name = P_NAME)
    params:
      sif = config['cactus_sif']
    log: "logs/cactus/hal_check.log"
    shell:
      """
      apptainer exec --bind $(pwd),{s_bind_paths} {params.sif} halStats {input} > {output}
      """
