"""
snakemake --rerun-triggers mtime -n final_figures
"""

rule final_figures:
    input:
      main_f2 = "../results/img/final/f_dem.pdf",
      main_f4 = "../results/img/final/f_whg.pdf",
      si_sf6 = "../results/img/final/sf_dem.pdf"

# requires imported results from RAD analysis
rule fig_dem:
    input:
      boots = "../results/demography/all_models_raw_boots.tsv",
      ci = "../results/demography/all_models_ci_boots.tsv",
      estimates = "../results/demography/all_models_estimates.tsv"
    output: "../results/img/final/f_dem.pdf"
    conda: "r_tidy"
    shell:
      """
      Rscript --vanilla R/figures/plot_main_demography.R
      """

rule fig_whg:
    input:
      p_pi =  "../results/img/R/p_pi.Rds",
      p_het = "../results/img/R/p_het_ind_bp_pheno.Rds",
      p_froh = "../results/img/R/p_f_rho_callable_ind.Rds",
      p_froh_cum_thresholds = "../results/img/R/p_cum_f_rho_callable_thresholds.Rds",
      p_froh_cum = "../results/img/R/p_cum_f_rho_callable_pheno.Rds",
      p_load_type = "../results/img/R/p_load_by_type_b.Rds",
      p_load_ind = "../results/img/R/p_load_by_ind_b.Rds"
    output: "../results/img/final/f_whg.pdf"
    conda: "r_tidy"
    shell:
      """
      Rscript --vanilla R/figures/plot_main_whg.R
      """

# requires imported results from RAD analysis
rule sup_fig_dem:
    input:
      boots = "../results/demography/all_models_raw_boots.tsv",
      ci = "../results/demography/all_models_ci_boots.tsv",
      estimates = "../results/demography/all_models_estimates.tsv"
    output: "../results/img/final/sf_dem.pdf"
    conda: "r_tidy"
    shell:
      """
      Rscript --vanilla R/figures/plot_si_demography.R
      """