"""
*auxillary rules for converions*

snakemake -n ../results/genotyping/filtered/mirang_filtered.ped
"""

rule vcf_to_plink:
    input: 
      vcf = "../results/genotyping/filtered/{file_base}.vcf.gz"
    output:
      pl_bed = "../results/genotyping/filtered/{file_base}.bed",
      pl_bim = "../results/genotyping/filtered/{file_base}.bim",
      pl_fam = "../results/genotyping/filtered/{file_base}.fam",
      pl_map = "../results/genotyping/filtered/{file_base}.map",
      pl_nosex = "../results/genotyping/filtered/{file_base}.nosex",
      pl_ped = "../results/genotyping/filtered/{file_base}.ped"
    params:
      geno_base = "../results/genotyping/filtered/"
    benchmark:
      "benchmark/conversion/vcf_to_plink_{file_base}.tsv"
    resources:
      mem_mb=15360
    container: c_popgen
    shell:
        """
        vcftools \
          --gzvcf {input.vcf} \
          --plink \
          --out {params.geno_base}{wildcards.file_base}

        plink \
          --file {params.geno_base}{wildcards.file_base} \
          --make-bed \
          --allow-extra-chr \
          --out {params.geno_base}{wildcards.file_base}
        """
