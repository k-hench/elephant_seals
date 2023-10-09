"""
snakemake --rerun-triggers mtime -n final_figures

snakemake --rulegraph final_figures | dot -Tsvg > ../results/img/control/dag_final_figures.svg
"""

rule final_figures:
    input:
      main_f2 = "../results/img/final/f_dem.pdf",
      main_f4 = "../results/img/final/f_whg.pdf",
      si_sf6 = "../results/img/final/sf_dem.pdf",
      si_sfX = "../results/img/final/sf_roh_length.pdf",
      si_sfxx = "../results/img/final/sf_het.pdf"

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

rule sup_fig_roh_length:
    input:
      p_froh_cum_thresholds = "../results/img/R/p_cum_f_rho_callable_thresholds.Rds"
    output: "../results/img/final/sf_roh_length.pdf"
    conda: "r_tidy"
    shell:
      """
      Rscript --vanilla R/figures/plot_si_roh_length.R
      """


rule sup_fig_heterozygosity:
    input:
      fai = "../data/genomes/filtered/mirang_filt.fa.gz.fai",
      spec_pops = expand( "../results/pop/inds_{spec}.pop", spec = BOTH_SPECS ),
      phenotypes = "../results/pop/group_pheno_labeled.pop",
      heterozygosity = "../results/het/win_het_ind_all_w1Mb_s250kb.tsv.gz",
      mirang_sex_chr = "../results/genomes/sex_chrom/mirang_sex_chrom.bed",
      sample_id = "../data/file_info.tsv",
      roh_data = expand( "../results/roh/bcftools/bed/max_callable/roh_max_{sample}_on_mirang.bed", sample = SAMPLES )
    output: "../results/img/final/sf_het.pdf"
    conda: "r_tidy"
    shell:
      """
      Rscript --vanilla R/figures/plot_si_roh_heterozygosity.R
      """