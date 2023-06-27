"""
*auxillary rules for converions*

snakemake --rerun-triggers mtime  -n ../results/genotyping/filtered/mirang_filtered.ped

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

rule create_species_pop:
    output: 
      pp = "../results/inds_{spec}.pop"
    params: 
      inds = lambda wc: SAMPLES_SPEC[wc.spec]
    shell:
      """
      echo "{params.inds}" | sed "s/ /\\n/g; s/'//g; s/\[//g; s/\]//g; /^[[:space:]]*$/d" > {output.pp}
      """

rule vcf_subset_species:
    input:
      vcf = "../results/genotyping/filtered/{file_base}.vcf.gz",
      inds = "../results/inds_{spec}.pop"
    output:
      vcf = "../results/genotyping/filtered/{file_base}_{spec}.vcf.gz"
    log:
      "logs/conversion/subset_{file_base}_{spec}.log"
    resources:
      mem_mb=15360
    container: c_popgen
    shell:
      """
      vcftools \
          --gzvcf {input.vcf} \
          --indv {input.inds} \
          --mac 1 \
          --recode \
          --stdout | \
          bgzip > {output.vcf} &> {log}
      
      tabix -p vcf {output.vcf} &>> {log}
      """