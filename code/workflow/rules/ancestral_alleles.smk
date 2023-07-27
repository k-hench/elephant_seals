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

rule merge_ancestral_hals:
    input:
      hal_anc = "../results/ancestral_allele/Anc52.hal",
      hal_mir = "../results/ancestral_allele/Anc56.hal"
    output:
      hal = "../results/ancestral_allele/mir_ancestral.hal",
      check = "../results/ancestral_allele/mir_ancestral_summary.txt"
    resources:
      mem_mb=15360
    container: c_cactus
    shell:
      """
      cp {input.hal_anc} {output.hal}
      halAppendSubtree {output.hal} {input.hal_mir} Anc56 Anc56 --merge --hdf5InMemory --hdf5InMemory
      halStats {output.hal} > {output.check}
      """

rule anc_allele_tsv:
    input:
      hal = "../results/ancestral_allele/mir_ancestral.hal"
    output:
      tsv = "../results/ancestral_allele/mirang_anc52_snps.tsv"
    resources:
      mem_mb=15360
    container: c_cactus
    shell:
      """
      halSnps {input.hal} mirang Anc52 --tsv {output.tsv}
      """

rule anc_tree:
    input:
      hal = "../results/ancestral_allele/mir_ancestral.hal"
    output:
      nwk = "../results/ancestral_allele/mirang_anc_tree.nwk"
    container: c_cactus
    shell:
      """
      halStats {input.hal} --tree  > {output.nwk}
      """