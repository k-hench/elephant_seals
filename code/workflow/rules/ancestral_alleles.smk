"""
snakemake --rerun-triggers mtime -n all_anc_allele
snakemake --rerun-triggers mtime   --use-singularity --singularity-args "--bind $CDATA" --use-conda -c 1 all_anc_allele

snakemake --jobs 100 \
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
  --jn job_aa.{name}.{jobid}.sh \
  -R all_anc_allele
"""

rule all_anc_allele:
    input: ""

rule extract_ancestral_hals:
    input:
      hal_anc = "../results/cactus/scratch/carnivora_set/tmp/steps-output/Anc52.hal",
      hal_mir = "../results/cactus/scratch/carnivora_set/tmp/steps-output/Anc56.hal"
    output:
      hal_anc = "../results/ancestral_allele/Anc52.hal",
      hal_mir = "../results/ancestral_allele/Anc56.hal"
    shell:
      """
      cp {input.hal_anc} {output.hal_anc}
      cp {input.hal_mir} {output.hal_mir}
      """