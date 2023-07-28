"""
snakemake --rerun-triggers mtime -n all_diversity
"""

rule all_diversity:
    input: ""

rule he_by_ind:
    input:
      vcf = "../results/genotyping/autosome/mirang_filtered_mirang_autosome.vcf.gz",
      inds = "../results/pop/inds_{spec}.pop"
    output:
      het = "../results/het/het_{spec}.tsv"
    container: c_popgen
    shell:
      """
      vcftools \
        --gzvcf {input.vcf} \
        --keep {input.inds} \
        --mac 1 \
        --het \
        --stdout > {output.het}
      """

wildcard_constraints:
    part="[^_]*"

rule all_bp_geno_by_spec:
    input:
      vcf = "../results/genotyping/filtered/partitions/mirang_all_bp_{part}_filtered.vcf.gz",
      inds = "../results/pop/inds_{spec}.pop"
    output:
      geno = "../results/genotyping/geno/partitions/{spec}_all_bp_{part}.geno"
    shell:
      """
      parseVCF -i minimal.vcf  | gzip > minimal.geno.gz
      """
    