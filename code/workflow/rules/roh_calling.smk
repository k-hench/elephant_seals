"""
snakemake --rerun-triggers mtime -n -R call_roh
# >>> needs to be run on ALL BP
"""

wildcard_constraints:
    het = "[0-9]*",
    whet = "[0-9]*",
    nsnp = "[0-9]*",
    wnsnp = "[0-9]*",
    leng = "[0-9]*",
    gap = "[0-9]*"

rule call_roh:
    input:
      expand( "../results/roh/bcftools/snp_based/bed/max_certain/roh_cert_{sample}_on_{ref}.bed", ref = GATK_REF[0], sample = SAMPLES ),
      plink_roh = expand( "../results/roh/plink/{file_base}_h{het}_wh{whet}_n{nsnp}_wn{wnsnp}_l{leng}_g{gap}", file_base = "mirang_filtered_all", het = [0], whet = [2], nsnp = [10], wnsnp = [50], leng = [10], gap = [1] )

'''
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
      SAMPLES=$(echo {params.samples} | sed 's/ /,/g')

      bcftools \
          roh {input.vcf} \
          -e - \
          -s $SAMPLES \
          -O rz > {output.roh} 2> {log}
    
      echo -e "-----------------" >> {log} 
      bcftools \
          roh {input.vcf} \
          -e - \
          -s $SAMPLES \
          -O sz > {output.roh_snps} 2>> {log}
      """

rule roh_to_bed:
    input:
      roh = expand("../results/roh/bcftools/{{ref}}_{part}_roh.tsv.gz", part =  GENOME_PARTITIONS )
    output:
      bed = "../results/roh/bcftools/bed/max_callable/roh_max_{sample}_on_{ref}.bed"
    shell:
      """
      for k in {input.roh}; do
        zcat $k | awk '$2=="{wildcards.sample}"{{print $3"\t"$4-1"\t"$5}}' >> {output.bed};
      done
      """

rule roh_max_certain:
    input:
      roh = "../results/roh/bcftools/bed/max_callable/roh_max_{sample}_on_{ref}.bed",
      cov_mask = "../results/qc/coverage/masks/{sample}_on_{ref}_binary_covmask.bed.gz"
    output:
      roh = "../results/roh/bcftools/bed/max_certain/roh_cert_{sample}_on_{ref}.bed",
    benchmark:
      'benchmark/roh/max_certain_roh_{sample}_on_{ref}.tsv'
    params:
      min_cov = 2
    resources:
      mem_mb=40960
    container: c_popgen
    shell:
      """
      zcat {input.cov_mask} |
        intersectBed -a stdin -b {input.roh} > {output.roh}
      """
'''

rule roh_plink:
    input:
      pl_bed = "../results/genotyping/plink/{file_base}.bed",
      pl_bim = "../results/genotyping/plink/{file_base}.bim",
      pl_fam = "../results/genotyping/plink/{file_base}.fam",
      pl_map = "../results/genotyping/plink/{file_base}.map",
      pl_nosex = "../results/genotyping/plink/{file_base}.nosex",
      pl_ped = "../results/genotyping/plink/{file_base}.ped"
    output:
      plink_dir = directory( "../results/roh/plink/{file_base}_h{het}_wh{whet}_n{nsnp}_wn{wnsnp}_l{leng}_g{gap}" )
    params:
      pl_base = "../results/genotyping/plink/{file_base}",
      out_dir = "../results/roh/plink",
      out_prefix = "{file_base}_h{het}_wh{whet}_n{nsnp}_wn{wnsnp}_l{leng}_g{gap}"
      #            "../results/roh/plink/{file_base}_h0_wh2_n10_wn50_l10_g1"
    resources:
      mem_mb=15360
    container: c_popgen
    shell:
      """
      mkdir -p {params.out_dir}/{params.out_prefix}
      plink --bfile {params.pl_base} \
        --out {params.out_dir}/{params.out_prefix}/{params.out_prefix} \
        --homozyg \
        --homozyg-window-snp {wildcards.wnsnp} \
        --homozyg-snp {wildcards.nsnp} \
        --homozyg-kb {wildcards.leng} \
        --homozyg-gap {wildcards.gap} \
        --homozyg-density 1 \
        --homozyg-window-missing 1 \
        --homozyg-het {wildcards.whet} \
        --homozyg-window-het {wildcards.het}
      """

rule roh_calling_bcftools_snps_only:
    input:
      vcf = "../results/genotyping/filtered/{ref}_filtered_{spec}.vcf.gz"
    output:
      roh = "../results/roh/bcftools/snp_based/{spec}_on_{ref}_roh.tsv.gz",
      roh_snps  = "../results/roh/bcftools/snp_based/{spec}_on_{ref}_roh_snps.tsv.gz"
    params:
      samples = lambda wc: expand( "{smp}", smp = SAMPLES_SPEC[wc.spec] )
    benchmark:
      "benchmark/roh/bcftools_{spec}_on_{ref}_snp.tsv"
    log:
      "logs/roh/bcftools_{spec}_on_{ref}_snp.log"
    resources:
      mem_mb=35840
    container: c_popgen
    shell:
      """
      SAMPLES=$(echo {params.samples} | sed 's/ /,/g')

      bcftools \
          roh {input.vcf} \
          -e - \
          -s $SAMPLES \
          -O rz > {output.roh} 2> {log}
    
      echo -e "-----------------" >> {log} 
      bcftools \
          roh {input.vcf} \
          -e - \
          -s $SAMPLES \
          -O sz > {output.roh_snps} 2>> {log}
      """

def get_spec_from_sample(smpl):
    return( seq_file_data["spec"][seq_file_data["sample_id"] == smpl].values[0] )

rule roh_to_bed_snps_only:
    input:
      roh = lambda wc: expand( "../results/roh/bcftools/snp_based/{spec}_on_{{ref}}_roh.tsv.gz", spec = get_spec_from_sample(wc.sample) )
    output:
      bed = "../results/roh/bcftools/snp_based/bed/max_callable/roh_max_{sample}_on_{ref}.bed"
    shell:
      """
      zcat {input.roh} | awk '$2=="{wildcards.sample}"{{print $3"\t"$4-1"\t"$5}}' >> {output.bed};
      """

rule roh_max_certain_snps_only:
    input:
      roh = "../results/roh/bcftools/snp_based/bed/max_callable/roh_max_{sample}_on_{ref}.bed",
      cov_mask = "../results/qc/coverage/masks/{sample}_on_{ref}_binary_covmask.bed.gz"
    output:
      roh = "../results/roh/bcftools/snp_based/bed/max_certain/roh_cert_{sample}_on_{ref}.bed",
    benchmark:
      'benchmark/roh/max_certain_roh_snps_{sample}_on_{ref}.tsv'
    params:
      min_cov = 2
    resources:
      mem_mb=40960
    container: c_popgen
    shell:
      """
      zcat {input.cov_mask} |
        intersectBed -a stdin -b {input.roh} > {output.roh}
      """