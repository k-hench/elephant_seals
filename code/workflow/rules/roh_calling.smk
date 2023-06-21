"""
snakemake --rerun-triggers mtime  -n ../results/roh/bcftools/mirang_02_roh.tsv.gz
# >>> needs to be run on ALL BP
"""

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