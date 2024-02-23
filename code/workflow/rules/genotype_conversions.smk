"""
*auxillary rules for converions*

snakemake --rerun-triggers mtime  -n ../results/genotyping/plink/mirang_filtered-mac2.ped

# python dummy wildcards
wc_dummy = type('MyObject', (object,), {})
wildcards = wc_dummy()
wildcards.spec = "all"
wildcards.spec = "mirang"
"""

# get all sample ids for an individual species
def get_samples_spec(spec):
    return( seq_file_data["sample_id"][seq_file_data["spec"] == spec].values )

def get_samples_spec_wc(wildcards):
    if( wildcards.spec == "all" ):
      return( seq_file_data["sample_id"].values )
    else:
      return( seq_file_data["sample_id"][seq_file_data["spec"] == wildcards.spec].values )

SAMPLES_ANG = sorted(list(set(get_samples_spec("mirang"))))
SAMPLES_LEO = sorted(list(set(get_samples_spec("mirleo"))))
SAMPLES_SPEC = {"mirang": SAMPLES_ANG, "mirleo": SAMPLES_LEO, "all": SAMPLES}

rule vcf_to_plink:
    input: 
      vcf = "../results/genotyping/filtered/{file_base}.vcf.gz"
    output:
      pl_bed = "../results/genotyping/plink/{file_base}.bed",
      pl_bim = "../results/genotyping/plink/{file_base}.bim",
      pl_fam = "../results/genotyping/plink/{file_base}.fam",
      pl_map = "../results/genotyping/plink/{file_base}.map",
      pl_nosex = "../results/genotyping/plink/{file_base}.nosex",
      pl_ped = "../results/genotyping/plink/{file_base}.ped"
    params:
      geno_base = "../results/genotyping/plink/"
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

rule create_species_pop:
    output: 
      pp = "../results/pop/inds_{spec}.pop"
    params: 
      inds = lambda wc: SAMPLES_SPEC[wc.spec]
    shell:
      """
      echo "{params.inds}" | sed "s/ /\\n/g; s/'//g; s/\[//g; s/\]//g; /^[[:space:]]*$/d" > {output.pp}
      """

rule vcf_subset_species:
    input:
      vcf = "../results/genotyping/filtered/{file_base}.vcf.gz",
      inds = "../results/pop/inds_{spec}.pop"
    output:
      vcf = "../results/genotyping/filtered/{file_base}_{spec}.vcf.gz"
    log:
      "logs/conversion/subset_{file_base}_{spec}.log"
    resources:
      mem_mb=15360
    container: c_popgen
    shell:
      """
      bcftools view \
        --samples-file {input.inds} \
        -Ov {input.vcf} \
        --min-ac 1:nonmajor | \
          bgzip > {output.vcf} 2> {log}
      
      tabix -p vcf {output.vcf} &>> {log}
      """

rule identify_chrX_scaffolds:
    input:
      size_m = "../results/genomes/mirang.size",
      size_z = "../results/genomes/zalcal.size",
      psl = "../results/psl/slim_mirang_on_zalcal.psl.gz"
    output:
      bed = "../results/genomes/sex_chrom/mirang_sex_chrom.bed"
    conda: "r_tidy"
    shell:
      """
      Rscript R/identify_scaffolds_on_x.R 
      """

rule vcf_filter_autosome:
    input:
      vcf = "../results/genotyping/filtered/{file_base}.vcf.gz",
      bed = "../results/genomes/sex_chrom/mirang_sex_chrom.bed"
    output:
      vcf = "../results/genotyping/autosome/{file_base}_autosome.vcf.gz"
    log:
      "logs/conversion/subset_{file_base}_autosome.log"
    resources:
      mem_mb=15360
    container: c_popgen
    shell:
      """
      vcftools \
          --gzvcf {input.vcf} \
          --exclude-bed {input.bed} \
          --recode \
          --stdout | \
          bgzip > {output.vcf} 2> {log}
      
      tabix -p vcf {output.vcf} &>> {log}
      """

rule vcf_filter_mac2:
    input:
      vcf = "../results/genotyping/filtered/{file_base}.vcf.gz",
    output:
      vcf = "../results/genotyping/filtered/{file_base}-mac2.vcf.gz"
    log:
      "logs/conversion/subset_{file_base}_mac2.log"
    resources:
      mem_mb=15360
    container: c_popgen
    shell:
      """
      vcftools \
          --gzvcf {input.vcf} \
          --mac 2 \
          --recode \
          --stdout | \
          bgzip > {output.vcf} 2> {log}
      
      tabix -p vcf {output.vcf} &>> {log}
      """
