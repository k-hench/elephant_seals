"""
# first create `carnivora_genome_and_timetree.tsv`
# to avoid errors raised by the function `get_accession`
snakemake -c 1 --use-conda ../data/carnivora_genome_and_timetree.tsv

snakemake -n -R ncbi_download

snakemake --dag -R  ncbi_download | dot -Tsvg > ../results/img/control/dag_download.svg

# interactive (online)
snakemake \
  --jobs 10 \
  --cores 10 \
  --latency-wait 30 \
  --use-conda \
  -R only_download

# batch-job submission (offline)
snakemake --jobs 70 \
  --latency-wait 30 \
  --use-conda \
  -p \
  --default-resources mem_mb=25600 \
  --cluster '
    qsub \
      -V -cwd \
      -P fair_share \
      -l idle=1 \
      -l si_flag=1 \
      -pe multislot {threads} \
      -l vf={resources.mem_mb}' \
      --jn job.{name}.{jobid}.sh \
      -R ncbi_download && mv job.* logs/
"""

WITH_CHICKEN = deepcopy( SPEC_ALL )
WITH_CHICKEN.insert( 0, "galgal" )

rule ncbi_download:
    message:
      """
      After the inital download from NCBI,
      the `arcgaz` genome needs to be replaced by
      the new dovetail assembly manually.
      """
    input: 
      '../results/img/qc/genomes_n50.svg',
      expand("../results/masking/{spec}_mask_check.tsv", spec = WITH_CHICKEN)

rule only_download:
    message:
      """
      bundeling of jobs that need online access
      (which might be limmited on cluster batch-job nodes)
      """
    input: 
      expand("../results/genomes/{spec}/{spec}.zip", spec = WITH_CHICKEN),
      expand("../results/genome_stats/{spec}.tsv", spec = WITH_CHICKEN)

checkpoint species_list:
    output: 
      species_list = "../data/carnivora_genome_and_timetree.tsv",
      short_label_tree = "../data/carnivora_short_labels.nwk",
      list_incl_chicken = "../data/carnivora_and_chicken.tsv"
    log:
      "../code/logs/r_species_list.log"
    container: None
    conda: "r_tidy"
    shell: 'Rscript R/compile_species_list.R 2> {log} 1> {log}'

def get_accession(wildcards, what):
    accessions = pd.read_table('../data/carnivora_and_chicken.tsv').set_index("spec", drop = False)
    if what == 'name':
        return accessions.loc[wildcards.spec, 'organism_name']
    elif what == 'accession':
        return accessions.loc[wildcards.spec, 'assembly_accession']
    elif what == 'repo':
        return accessions.loc[wildcards.spec, 'repo']

rule genome_stats:
    output: '../results/genome_stats/{spec}.tsv'
    params:
      name = lambda wc: get_accession(wc, what = "name"),
      accnr = lambda wc: get_accession(wc, what = "accession"),
      repo = lambda wc: get_accession(wc, what = "repo")
    conda:
      'ncbi_datasets'
    container: None
    shell:
      """
      datasets summary \
        genome accession "{params.accnr}" \
        --assembly-source {params.repo}  \
        --as-json-lines | \
        dataformat tsv \
        genome --fields \
        organism-name,accession,assminfo-name,annotinfo-name,assmstats-scaffold-n50,assmstats-contig-l50,assmstats-total-sequence-len,annotinfo-release-date \
        > {output}
      """

rule download_genome:
    output: '../results/genomes/{spec}/{spec}.zip'
    params:
      name = lambda wc: get_accession(wc, what = "name"),
      accnr = lambda wc: get_accession(wc, what = "accession")
    conda:
      'ncbi_datasets'
    container: None
    log:
      "logs/genome_dl/genome_dl_{spec}.log"
    shell:
      """
      mkdir -p ../results/genomes/{wildcards.spec}

      datasets download \
        genome \
        accession {params.accnr} \
        --filename {output} \
        --reference \
        --include genome \
         2> {log} 1> {log}
      """

rule repack_genome:
    input: 
      zp = '../results/genomes/{spec}/{spec,[a-z]+}.zip'
    output:
      gz = '../data/genomes/{spec,[a-z]+}.fa.gz'
    params:
      name = lambda wc: get_accession(wc, what = "name"),
      accnr = lambda wc: get_accession(wc, what = "accession")
    container: None
    conda:
      "map_align"
    log:
      "logs/genome_uz/genome_uz_{spec}.log"
    shell:
      """
      unzip -p \
        {input.zp} \
        ncbi_dataset/data/{params.accnr}/{params.accnr}*.fna | \
        bgzip > {output.gz}
      """

rule check_if_masked:
    input: '../data/genomes/{spec}.fa.gz'
    output: '../results/masking/{spec}_mask_check.tsv'
    log:
      "logs/mask_check/mask_check_{spec}.log"
    container: None
    conda:
      "map_align"
    shell:
      """
      mkdir -p ../results/masking/
      
      # check if grep fails (no match)
      zgrep -v "^>" {input} | grep -q '[atgc]' || GREPERR=$? 
      echo "error code from grep: "$GREPERR

      if [ $GREPERR -eq 1 ]; then 
        echo -e "{wildcards.spec}\t0\tunmasked" > {output}
      else 
        echo -e "{wildcards.spec}\t1\tmasked" > {output}
      fi
      """

rule stat_plots:
    input: 
      stats = expand("../results/genome_stats/{spec}.tsv", spec = SPEC_ALL),
      genomes = expand("../data/genomes/{spec}.fa.gz", spec = SPEC_ALL)
    output: "../results/img/qc/genomes_n50.svg"
    container: None
    conda: "r_tidy"
    log:
      "logs/r_genome_stats.log"
    shell:
      """
      Rscript R/genome_stats.R 2> {log} 1> {log}
      """