"""
snakemake --rerun-triggers mtime  -n -R call_roh
# >>> needs to be run on ALL BP
"""

rule call_roh:
    input: expand( "../results/roh/bcftools/{ref}_{part}_roh_snps.tsv.gz", ref = GATK_REF[0], part =  GENOME_PARTITIONS )

rule roh_calling_bcftools:
    input: 
      vcf = "../results/genotyping/filtered/{ref}_all_bp_{part}_filtered.vcf.gz"
    output:
      roh = "../results/roh/bcftools/{ref}_{part}_roh.tsv.gz",
      roh_snps  = "../results/roh/bcftools/{ref}_{part}_roh_snps.tsv.gz"
    params:
      samples = expand( "{smp}", smp = SAMPLES )
    benchmark:
      "benchmark/roh/bcftools_{ref}_{part}.tsv"
    resources:
      mem_mb=15360
    container: c_popgen
    shell:
      """
      bcftools \
          roh {input.vcf} \
          -e - \
          -s {params.samples} \
          -O rz > {output.roh}
    
      bcftools \
          roh {input.vcf} \
          -e - \
          -s {params.samples} \
          -O sz > {output.roh_snps}
      """