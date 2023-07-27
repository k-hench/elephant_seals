"""
snakemake --rerun-triggers mtime  -n all_demography
very much based upon the workshop "speciation genomics"
by Joana Meier and Mark Ravinet
https://speciationgenomics.github.io/fastsimcoal2/

snakemake --jobs 100 \
  --latency-wait 30 \
  -p \
  --default-resources mem_mb=51200 threads=1 \
  --use-singularity \
  --singularity-args "--bind $CDATA" \
  --use-conda \
  --rerun-triggers mtime \
  --cluster '
    qsub \
      -V -cwd \
      -P fair_share \
      -l idle=1 \
      -l si_flag=1 \
      -pe multislot {threads} \
      -l vf={resources.mem_mb}' \
  --jn job_fs.{name}.{jobid}.sh \
  -R all_demography
"""

DEM_TYPES = [ "bot06-lgm", "bot06-nes", "bot10-lgm", "bot10-nes", "null-lgm", "null-nes" ]
DEM_N = [ str(x + 1).zfill(3) for x in np.arange(100) ]
BOOTSTRAP_N = [ str(x + 1).zfill(2) for x in np.arange(50) ]

rule all_demography:
    input: 
      sfs_prev = expand( "../results/demography/preview/prev_{spec}_on_{ref}.txt" , ref = "mirang", spec = "mirang" ),
      sfs_dir = expand( "../results/demography/sfs/{spec}_on_{ref}" , ref = "mirang", spec = "mirang" ),
      fs_iter = expand( "../results/demography/fastsimcoal/{spec}_on_{ref}/{fs_run}/bestrun/{spec}_on_{ref}_{fs_run}.lhoods", ref = "mirang", spec = "mirang", fs_run = DEM_TYPES ),
      bs_idx = expand( "../results/demography/bootstrap/{spec}_on_{ref}_bs_{idx}", ref = "mirang", spec = "mirang", idx = BOOTSTRAP_N ),
      bs_best = expand( "../results/demography/fastsimcoal/{spec}_on_{ref}/{fs_run}/bootstrap/bs_{idx}/all_lhoods.tsv", ref = "mirang", spec = "mirang", fs_run = DEM_TYPES, idx = BOOTSTRAP_N )

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
      vcf = "../results/genotyping/autosome/{ref}_filtered_{spec}_autosome.vcf.gz"
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
      vcf = "../results/genotyping/autosome/{ref}_filtered_{spec}_autosome.vcf.gz"
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
    container: c_sim
    shell:
      """
      cp {params.obs} {output.obs}
      cd {params.basedir}/bestrun

      # Run fastsimcoal 20 times (in reality better 100 times) to get the likelihood of the observed SFS under the best parameter values with 1 mio simulated SFS.
      for i in {{1..100}}; do
        echo $i >> loop.log
        fsc27093 -i {params.prefix}_maxL.par -n 1000000 -m -q -0
        # Fastsimcoal will generate a new folder called {params.prefix}_maxL and write files in there

        # collect the lhood values
        sed -n '2,3p' {params.prefix}_maxL/{params.prefix}_maxL.lhoods  >> {params.prefix}.lhoods

        # delete the folder with results
        rm -r {params.prefix}_maxL/
      done
      """

rule bootstrap_prep:
    input:
      vcf = "../results/genotyping/filtered/{ref}_filtered_{spec}.vcf.gz"
    output:
      head = "../results/genotyping/bootstrap/{ref}_filtered_{spec}.header",
      all_sites = temp( "../results/genotyping/bootstrap/{ref}_filtered_{spec}.allSites" )
    params:
      n_sites = 11933,
      block_base = "../results/genotyping/bootstrap/{ref}_filtered_{spec}.sites."
    shell:
      """
      zgrep -v "^#" {input.vcf} > {output.all_sites}

      # Get the header
      zgrep "^#" {input.vcf} > {output.head}

      # get 100 files with {params.n_sites} sites each
      split -l {params.n_sites} {output.all_sites} {params.block_base}
      """

rule bootstrap_vcf:
    input:
      head =  "../results/genotyping/bootstrap/{ref}_filtered_{spec}.header"
    output:
      vcf_bs = "../results/genotyping/bootstrap/{ref}_filtered_{spec}_bs_{idx}.vcf.gz"
    params:
      block_base = "../results/genotyping/bootstrap/{ref}_filtered_{spec}.sites.",
      vcf_base = "../results/genotyping/bootstrap/{ref}_filtered_{spec}_bs_{idx}.vcf"
    container: c_popgen
    shell:
      """
      cp {input.head} {params.vcf_base}
      
      # Randomly add 100 blocks
        for r in {{1..100}}; do
          cat `shuf -n1 -e {params.block_base}*` >> {params.vcf_base}
        done
      
      # Compress the vcf file again
      bgzip {params.vcf_base}
      """

rule bootstrap_sfs:
    input:
      pp = "../results/pop/inds_{spec}.pop2",
      vcf_bs = "../results/genotyping/bootstrap/{ref}_filtered_{spec}_bs_{idx}.vcf.gz"
    output:
      sfs_dir = directory( "../results/demography/bootstrap/{spec}_on_{ref}_bs_{idx}" )
    params:
      n_haplo = 40
    container: c_sim
    shell:
      """
      # Make an SFS from the new bootstrapped file
      easySFS.py -i {input.vcf_bs} -p {input.pp} -a -f --proj {params.n_haplo} -o {output.sfs_dir}
      """

rule bootstrap_fastsimcoal:
    input:
      sfs_dir = "../results/demography/bootstrap/{spec}_on_{ref}_bs_{idx}",
      tpl = "../data/templates/tpl/{fs_run}.tpl",
      est =  "../data/templates/est/{fs_run}.est"
    output:
      fs_dir = temp( directory( "../results/demography/fastsimcoal/{spec}_on_{ref}/{fs_run}/bootstrap/bs_{idx}/bs{idx}_{fs_run}_{iter}" ) )
    params:
      runs = "{fs_run}",
      prefix = "{spec}_on_{ref}_{fs_run}",
      basedir = "../results/demography/fastsimcoal/{spec}_on_{ref}/{fs_run}/bootstrap/bs_{idx}",
      obs = "../results/demography/bootstrap/{spec}_on_{ref}_bs_{idx}/fastsimcoal2/{spec}_MAFpop0.obs"
    resources:
      mem_mb=15360
    threads: 4
    container: c_sim
    shell:
      """
      mkdir -p {output.fs_dir}
      # this time SFS from bootsraped file
      cp {params.obs} {output.fs_dir}/{params.prefix}_MAFpop0.obs
      # double check that these remain the same
      cp {input.tpl} {output.fs_dir}/{params.prefix}.tpl
      cp {input.est} {output.fs_dir}/{params.prefix}.est
      cd {output.fs_dir}

      # n: number of simulations, m: minor sfs, M: max.likelihood, L: number of ECM cycles (loops)
      # q: quiet, w: tolerance for brent optimization, x: no arlequin output
      # C: min. observed SFS count, c: cores
      fsc27093 -t {params.prefix}.tpl -n 100000 -m -e {params.prefix}.est -M -L 40 -q -w 0.01 --foldedSFS -x -C 5 --nosingleton -c 4
      """

rule bootrap_best_run:
    input:
      all_runs = expand( "../results/demography/fastsimcoal/{{spec}}_on_{{ref}}/{{fs_run}}/bootstrap/bs_{{idx}}/bs{{idx}}_{{fs_run}}_{iter}", iter = DEM_N )
    output:
      all_lhoods = "../results/demography/fastsimcoal/{spec}_on_{ref}/{fs_run}/bootstrap/bs_{idx}/all_lhoods.tsv"
    params:
      runs = "{fs_run}",
      prefix = "{spec}_on_{ref}_{fs_run}",
      basedir = "../results/demography/fastsimcoal/{spec}_on_{ref}/{fs_run}/bootstrap/bs_{idx}"
    shell:
      """
      cd {params.basedir}
      FLS=$( ls bs{wildcards.idx}_{params.runs}_*/{params.prefix}/{params.prefix}.bestlhoods )

      echo -e "RUN\tMaxEstLhood\tMaxObsLhood\tDELTA_OBS_EST" > all_lhoods.tsv
      for k in $FLS; do
        RUNNR=$(echo $k | sed "s=/.*==g; s/bs{wildcards.idx}_{params.runs}_//")
        awk -v r="$RUNNR" 'NR==2{{print r"\t"$(NF-1)"\t"$NF"\t"$(NF-1)-$NF}}' $k >> all_lhoods.tsv
      done

      BEST_RUN=$(sort -k 4 all_lhoods.tsv  | head -n 1 | cut -f 1)

      cp -r bs{wildcards.idx}_{params.runs}_${{BEST_RUN}}/{params.prefix} ./bestrun
      """

wildcard_constraints:
    fs_run="[^_]*"