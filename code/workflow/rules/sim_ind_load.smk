"""
snakemake --rerun-triggers mtime -n all_ml_sim
snakemake --rerun-triggers mtime -c 10 --jobs 10 --use-singularity --singularity-args "--bind $CDATA" --use-conda  all_ml_sim

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
  -R all_ml_sim
"""

SIM_INDS = [ "i" + str(x) for x in np.arange(20)]
SIM_IDX = [ 1 + x for x in np.arange(100)]

rule all_ml_sim:
    input:
     s = expand( "../results/simulations/selection_coefficient/sim_{sim_idx}_s.tsv.gz", sim_idx = SIM_IDX ),
     tally = expand( "../results/simulations/tally/sim_{sim_idx}_load_by_ind.tsv", sim_idx = SIM_IDX )

rule sim_selection_coefficients:
    input: 
      vcf = "../results/simulations/NES_SLiM_VCFs/NES_sim_Sample_{sim_idx}.vcf"
    output:
      bed = "../results/simulations/selection_coefficient/sim_{sim_idx}_s.tsv.gz"
    resources:
      mem_mb=5120
    container: c_ml
    shell:
      """
      cat {input.vcf} | \
        SnpSift extractFields /dev/stdin CHROM POS S DOM | \
        gzip > {output.bed}
      """

rule sim_tally_load:
    input: 
      vcf = "../results/simulations/NES_SLiM_VCFs/NES_sim_Sample_{sim_idx}.vcf"
    output:
      tsv = "../results/simulations/tally/sim_{sim_idx}_load_by_ind.tsv"
    resources:
      mem_mb=5120
    container: c_ml
    shell:
      """
      S_IDX=1
      echo -e "sample\tn_masked_strong\tn_masked_intermediate\tn_masked_weak\tn_expressed_strong\tn_expressed_intermediate\tn_expressed_weak" > {output.tsv}
      for k in {{0..19}}; do
        ./sh/tally_load_ind.sh {input.vcf} i${{k}} >> {output.tsv}
      done 
      """
