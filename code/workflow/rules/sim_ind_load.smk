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
    input: expand( "../results/simulations/by_ind/{load_type}/sim_{sim_idx}/{sim_ind}_{load_type}.bed.gz", load_type = ["expressed", "masked"], sim_idx = SIM_IDX, sim_ind = SIM_INDS )


rule sim_masked_load:
    input: 
      vcf = "../results/simulations/NES_SLiM_VCFs/NES_sim_Sample_{sim_idx}.vcf"
    output:
      bed = "../results/simulations/by_ind/masked/sim_{sim_idx}/{sim_ind}_masked.bed.gz"
    resources:
      mem_mb=5120
    container: c_ml
    shell:
      """
      # heterozygous (masked load)
      zcat {input.vcf} | \
        SnpSift filter "( isHet(GEN[{wildcards.sim_ind}].GT) )" | \
        grep -v "^##" | \
        awk -v OFS="\t" -v s="{wildcards.sim_ind}" '{{if(NR==1){{ for (i=1; i<=NF; ++i) {{ if ($i ~ s) c=i }} }} {{print $1,$2,$2,$c}} }}' | \
        sed 's/POS\tPOS/FROM\tTO/' | \
        gzip > {output.bed}
      """

rule sim_expressed_load:
    input: 
      vcf = "../results/simulations/NES_SLiM_VCFs/NES_sim_Sample_{sim_idx}.vcf"
    output:
      bed = "../results/simulations/by_ind/expressed/sim_{sim_idx}/{sim_ind}_expressed.bed.gz"
    resources:
      mem_mb=5120
    container: c_ml
    shell:
      """
      # homozygous for derived allele (expressed load - ancestral def)
      zcat {input.vcf} | \
        SnpSift filter "((isVariant(GEN[{wildcards.sim_ind}].GT)) & (isHom(GEN[{wildcards.sim_ind}].GT)))" | \
        grep -v "^##" | \
        awk -v OFS="\t" -v s="{wildcards.sim_ind}" '{{if(NR==1){{ for (i=1; i<=NF; ++i) {{ if ($i ~ s) c=i }} }} {{print $1,$2,$2,$c}} }}' | \
        sed 's/POS\tPOS/FROM\tTO/' | \
        gzip > {output.bed}
      """