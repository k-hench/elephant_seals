"""
snakemake --rerun-triggers mtime -n -R call_roh
# >>> needs to be run on ALL BP
"""

rule call_roh:
    input:
      expand( "../results/roh/bcftools/{ref}_{part}_roh_snps.tsv.gz", ref = GATK_REF[0], part =  GENOME_PARTITIONS ),
      expand( "../results/roh/bcftools/snp_based/{ref}_roh_snps.tsv.gz", ref = GATK_REF[0] )

rule roh_calling_bcftools:
    input:
      vcf = "../results/genotyping/filtered/partitions/{ref}_all_bp_{part}_filtered.vcf.gz"
    output:
      roh = "../results/roh/bcftools/{ref}_{part}_roh.tsv.gz",
      roh_snps  = "../results/roh/bcftools/{ref}_{part}_roh_snps.tsv.gz"
    params:
      samples = expand( "{smp}", smp = SAMPLES )
    benchmark:
      "benchmark/roh/bcftools_{ref}_{part}.tsv"
    log:
      "logs/roh/bcftools_{ref}_{part}.log"
    resources:
      mem_mb=35840
    container: c_popgen
    shell:
      """
      bcftools \
          roh {input.vcf} \
          -e - \
          -s {params.samples} \
          -O rz > {output.roh} 2> {log}
    
      echo -e "-----------------" >> {log} 
      bcftools \
          roh {input.vcf} \
          -e - \
          -s {params.samples} \
          -O sz > {output.roh_snps} 2>> {log}
      """

rule roh_calling_bcftools_snps_only:
    input:
      vcf = "../results/genotyping/filtered/partitions/{ref}_filtered.vcf.gz"
    output:
      roh = "../results/roh/snp_based/bcftools/{ref}_roh.tsv.gz",
      roh_snps  = "../results/roh/bcftools/snp_based/{ref}_roh_snps.tsv.gz"
    params:
      samples = expand( "{smp}", smp = SAMPLES )
    benchmark:
      "benchmark/roh/bcftools_{ref}_snp.tsv"
    log:
      "logs/roh/bcftools_{ref}_snp.log"
    resources:
      mem_mb=35840
    container: c_popgen
    shell:
      """
      bcftools \
          roh {input.vcf} \
          -e - \
          -s {params.samples} \
          -O rz > {output.roh} 2> {log}
    
      echo -e "-----------------" >> {log} 
      bcftools \
          roh {input.vcf} \
          -e - \
          -s {params.samples} \
          -O sz > {output.roh_snps} 2>> {log}
      """