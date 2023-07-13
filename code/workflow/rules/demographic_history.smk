"""
snakemake --rerun-triggers mtime  -n all_demography
"""

DEM_TYPES = [ "bot06-lgm", "bot06-nes", "bot10-lgm", "bot10-nes", "null-lgm", "null-nes" ]
DEM_N = [ str(x + 1).zfill(3) for x in np.arange(100) ]

rule all_demography:
    input: 
      sfs_prev = expand( "../results/demography/preview/prev_{spec}_on_{ref}.txt" , ref = "mirang", spec = "mirang" ),
      sfs_dir = expand( "../results/demography/sfs/{spec}_on_{ref}" , ref = "mirang", spec = "mirang" ),
      fs_iter = expand( "../results/demography/fastsimcoal/{spec}_on_{ref}/{fs_run}/bestrun/{spec}_on_{ref}_{fs_run}.AIC", ref = "mirang", spec = "mirang", fs_run = DEM_TYPES )

rule create_pop2_files:
    input:
      pp = "../results/pop/inds_{spec}.pop"
    output:
      pp = "../results/pop/inds_{spec}.pop2"
    shell:
      """
      awk -v pop={wildcards.spec} '{{print $1"\t"pop}}' {input.pp} > {output.pp}
      """

rule estimate_projections:
    input:
      pp = "../results/pop/inds_{spec}.pop2",
      vcf = "../results/genotyping/filtered/{ref}_filtered_{spec}.vcf.gz"
    output:
      prev = "../results/demography/preview/prev_{spec}_on_{ref}.txt"
    container: c_sim
    shell:
      """
      easySFS.py -i {input.vcf} -p {input.pp} -a -f --preview > {output.prev}
      """

rule create_sfs:
    input:
      pp = "../results/pop/inds_{spec}.pop2",
      vcf = "../results/genotyping/filtered/{ref}_filtered_{spec}.vcf.gz"
    output:
      sfs_dir = directory( "../results/demography/sfs/{spec}_on_{ref}/" )
    params:
      n_haplo = 40
    container: c_sim
    shell:
      """
      easySFS.py -i {input.vcf} -p {input.pp} -a -f --proj {params.n_haplo} -o {output.sfs_dir}
      """

rule run_fastsimcoal:
    input:
      sfs_dir = "../results/demography/sfs/{spec}_on_{ref}",
      tpl = "../data/templates/tpl/{fs_run}.tpl",
      est =  "../data/templates/est/{fs_run}.est"
    output:
      fs_dir = temp( directory( "../results/demography/fastsimcoal/{spec}_on_{ref}/{fs_run}/{fs_run}_{iter}" ) )
    params:
      obs = "../results/demography/sfs/{spec}_on_{ref}/fastsimcoal2/{spec}_MAFpop0.obs",
      prefix = "{spec}_on_{ref}_{fs_run}"
    resources:
      mem_mb=15360
    threads: 4
    container: c_sim
    shell:
      """
      mkdir -p {output.fs_dir}
      cp {params.obs} {output.fs_dir}/{params.prefix}_MAFpop0.obs
      cp {input.tpl} {output.fs_dir}/{params.prefix}.tpl
      cp {input.est} {output.fs_dir}/{params.prefix}.est
      cd {output.fs_dir}

      # n: number of simulations, m: minor sfs, M: max.likelihood, L: number of ECM cycles (loops)
      # q: quiet, w: tolerance for brent optimization, x: no arlequin output
      # C: min. observed SFS count, c: cores
      fsc27093 -t {params.prefix}.tpl -n 100000 -m -e {params.prefix}.est -M -L 40 -q -w 0.01 --foldedSFS -x -C 5 --nosingleton -c 4
      """

rule collect_best_fsc_run:
    input:
      all_runs = expand( "../results/demography/fastsimcoal/{{spec}}_on_{{ref}}/{{fs_run}}/{{fs_run}}_{iter}", iter = DEM_N )
    output:
      all_lhoods = "../results/demography/fastsimcoal/{spec}_on_{ref}/{fs_run}/all_lhoods.tsv"
    params:
      runs = "{fs_run}",
      prefix = "{spec}_on_{ref}_{fs_run}",
      basedir = "../results/demography/fastsimcoal/{spec}_on_{ref}/{fs_run}"
    shell:
      """
      cd {params.basedir}
      FLS=$( ls {params.runs}_*/{params.prefix}/{params.prefix}.bestlhoods )

      echo -e "RUN\tMaxEstLhood\tMaxObsLhood\tDELTA_OBS_EST" > all_lhoods.tsv
      for k in $FLS; do
        RUNNR=$(echo $k | sed "s=/.*==g; s/{params.runs}_//")
        awk -v r="$RUNNR" 'NR==2{{print r"\t"$(NF-1)"\t"$NF"\t"$(NF-1)-$NF}}' $k >> all_lhoods.tsv
      done

      BEST_RUN=$(sort -k 4 all_lhoods.tsv  | head -n 1 | cut -f 1)

      cp -r {params.runs}_${{BEST_RUN}}/{params.prefix} ./bestrun
      """

rule calculate_aic:
    input:
      all_lhoods = "../results/demography/fastsimcoal/{spec}_on_{ref}/{fs_run}/all_lhoods.tsv",
      tpl = "../data/templates/tpl/{fs_run}.tpl",
      est =  "../data/templates/est/{fs_run}.est"
    output:
      aic = "../results/demography/fastsimcoal/{spec}_on_{ref}/{fs_run}/bestrun/{spec}_on_{ref}_{fs_run}.AIC"
    params:
      runs = "{fs_run}",
      prefix = "{spec}_on_{ref}_{fs_run}",
      basedir = "../results/demography/fastsimcoal/{spec}_on_{ref}/{fs_run}"
    conda: "r_tidy"
    shell:
      """
      cp {input.tpl} {params.basedir}/bestrun/{params.prefix}.tpl
      cp {input.est} {params.basedir}/bestrun/{params.prefix}.est
      cd {params.basedir}/bestrun
      Rscript {base_dir}/code/R/calculateAIC_kh.R {params.prefix}
      """

rule likelihood_ditributions_bestrun:
    input:
      aic = "../results/demography/fastsimcoal/{spec}_on_{ref}/{fs_run}/bestrun/{spec}_on_{ref}_{fs_run}.AIC",
      sfs_dir = "../results/demography/sfs/{spec}_on_{ref}"
    output:
      obs = "../results/demography/fastsimcoal/{spec}_on_{ref}/{fs_run}/bestrun/{spec}_on_{ref}_{fs_run}_maxL_MAFpop0.obs",
      lhoods = "../results/demography/fastsimcoal/{spec}_on_{ref}/{fs_run}/bestrun/{spec}_on_{ref}_{fs_run}.lhoods",
      looplog = temp("../results/demography/fastsimcoal/{spec}_on_{ref}/{fs_run}/bestrun/loop.log" )
    params:
      runs = "{fs_run}",
      prefix = "{spec}_on_{ref}_{fs_run}",
      basedir = "../results/demography/fastsimcoal/{spec}_on_{ref}/{fs_run}",
      obs = "../results/demography/sfs/{spec}_on_{ref}/fastsimcoal2/{spec}_MAFpop0.obs"
    threads: 4
    container: c_sim
    shell:
      """
      cp {params.obs} {output.obs}
      cd {params.basedir}/bestrun

      # Run fastsimcoal 20 times (in reality better 100 times) to get the likelihood of the observed SFS under the best parameter values with 1 mio simulated SFS.
      for i in {1..100}; do
        echo $i >> loop.log
        fsc27093 -i {params.prefix}_maxL.par -n 1000000 -m -q -0
        # Fastsimcoal will generate a new folder called {params.prefix}_maxL and write files in there

        # collect the lhood values
        sed -n '2,3p' {params.prefix}_maxL/{params.prefix}_maxL.lhoods  >> {params.prefix}.lhoods

        # delete the folder with results
        rm -r {params.prefix}_maxL/
      done
      """