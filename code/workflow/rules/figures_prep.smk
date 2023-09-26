"""
snakemake --rerun-triggers mtime -n all_fig_prep
"""
BOTH_SPECS = ["mirang", "mirleo"]

rule all_fig_prep:
    input:
      dem_boots = "../results/demography/all_models_raw_boots.tsv",
      dem_ci = "../results/demography/all_models_ci_boots.tsv",
      dem_estimates = "../results/demography/all_models_estimates.tsv",
      p_pi =  "../results/img/R/p_pi.Rds",
      p_het = "../results/img/R/p_het_ind_bp_pheno.Rds",
      p_froh = "../results/img/R/p_f_rho_callable_ind.Rds",
      p_froh_cum_thresholds = "../results/img/R/p_cum_f_rho_callable_thresholds.Rds",
      p_froh_cum = "../results/img/R/p_cum_f_rho_callable_pheno.Rds",
      p_load_type = "../results/img/R/p_load_by_type_b.Rds",
      p_load_ind = "../results/img/R/p_load_by_ind_b.Rds"

# --- demography figures -------------
# requires results imported from RAD analysis
rule data_dem:
    input:
      whg = "../results/checkpoints/deomgraphy_done.check",
      rad_est = "../results/RAD/demography/FSC_parameterEstiamtes.csv",
      rad_boot = expand( "../results/RAD/demography/Non_Param_boots_{tag}Model.txt", tag = ["Bott10", "Bott6_Best", "Null"]),
      rad_boot_extra = "../results/RAD/demography/NonParam_900ExtraBoots.txt" 
    output:
      boots = "../results/demography/all_models_raw_boots.tsv",
      ci = "../results/demography/all_models_ci_boots.tsv",
      estimates = "../results/demography/all_models_estimates.tsv"
    conda: "r_tidy"
    shell:
      """
      Rscript --vanilla R/fastsimcoal_data_merge_whg_rad.R
      """

# --- whole genome figures -------------

rule plotprep_whg_pi:
    input: expand( "../results/pi/mirang_pi_dxy_{part}.tsv.gz", part = GENOME_PARTITIONS )
    output: "../results/img/R/p_pi.Rds"
    conda: "r_tidy"
    shell:
      """
      Rscript --vanilla R/figures/plot_main_whg_prep_pi.R
      """

rule plot_data_whg_het_ind:
    input:
      fai = "../data/genomes/filtered/mirang_filt.fa.gz.fai",
      g012 = expand( "../results/genotyping/012/mirang_filtered_{spec}_012.tsv.gz", spec = BOTH_SPECS),
      g012_pos =  expand( "../results/genotyping/012/mirang_filtered_{spec}_012.012.pos", spec = BOTH_SPECS ),
      g012_ind = expand( "../results/genotyping/012/mirang_filtered_{spec}_012.012.indv", spec = BOTH_SPECS )
    output: "../results/het/win_het_ind_all_w1Mb_s250kb.tsv.gz"
    conda: "r_tidy"
    shell:
      """
      Rscript --vanilla R/het_data_by_ind.R
      """

rule plot_data_whg_het_spec:
    input:
      fai = "../data/genomes/filtered/mirang_filt.fa.gz.fai",
      g012 = expand( "../results/genotyping/012/mirang_filtered_{spec}_012.tsv.gz", spec = BOTH_SPECS ),
      g012_pos =  expand( "../results/genotyping/012/mirang_filtered_{spec}_012.012.pos", spec = BOTH_SPECS )
    output: "../results/het/win_het_by_spec_w1Mb_s250kb.tsv.gz"
    conda: "r_tidy"
    shell:
      """
      Rscript --vanilla R/het_data_by_spec.R
      """

rule plotprep_whg_het:
    input:
      file_info = "../data/file_info.tsv",
      ind = "../results/het/win_het_ind_all_w1Mb_s250kb.tsv.gz",
      spec = "../results/het/win_het_by_spec_w1Mb_s250kb.tsv.gz"
    output:
      p_pheno = "../results/img/R/p_het_ind_bp_pheno.Rds"
    conda: "r_tidy"
    shell:
      """
      Rscript --vanilla R/figures/plot_main_whg_prep_het.R
      """

rule plotprep_whg_froh:
    input:
      fai = "../data/genomes/filtered/mirang_filt.fa.gz.fai",
      spec_pops = expand( "../results/pop/inds_{spec}.pop", spec = BOTH_SPECS ),
      roh_data = expand( "../results/roh/bcftools/bed/max_callable/roh_max_{sample}_on_mirang.bed", sample = SAMPLES ),
      mirang_sex_chr = "../results/genomes/sex_chrom/mirang_sex_chrom.bed",
      meverage_masks = expand( "../results/qc/coverage/masks/{sample}_on_mirang_binary_covmask.bed.gz", sample = SAMPLES ),
      file_info = "../data/file_info.tsv"
    output:
      p_froh = "../results/img/R/p_f_rho_callable_ind.Rds",
      p_froh_cum_thresholds = "../results/img/R/p_cum_f_rho_callable_thresholds.Rds",
      p_froh_cum = "../results/img/R/p_cum_f_rho_callable_pheno.Rds"
    conda: "r_tidy"
    shell:
      """
      Rscript --vanilla R/figures/plot_main_whg_prep_roh.R
      """

rule plotprep_whg_load:
    input:
      file_info = "../data/file_info.tsv",
      load_data = expand( "../results/mutation_load/snp_eff/by_ind/{load_type}/{sample}_{load_type}.bed.gz", load_type = ["masked", "expressed", "fixed"], sample = SAMPLES ),
      load_roh_data = expand( "../results/mutation_load/snp_eff/by_ind/{load_type}_in_roh/{sample}_{load_type}_in_roh.bed.gz", load_type = ["masked", "expressed", "fixed"], sample = SAMPLES ),
      load_data_anc = expand("../results/mutation_load/snp_eff/by_ind/{load_type}_anc/{sample}_{load_type}_anc.bed.gz", load_type = ["expressed", "fixed"], sample = SAMPLES ),
      load_data_anc_roh = expand("../results/mutation_load/snp_eff/by_ind/{load_type}_anc_in_roh/{sample}_{load_type}_anc_in_roh.bed.gz", load_type = ["expressed", "fixed"], sample = SAMPLES )
    output:
      p_load_type = "../results/img/R/p_load_by_type_b.Rds",
      p_load_ind = "../results/img/R/p_load_by_ind_b.Rds"
    conda: "r_tidy"
    shell:
      """
      Rscript --vanilla R/figures/plot_main_whg_prep_roh.R
      """
