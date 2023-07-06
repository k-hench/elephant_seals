"""
snakemake --rerun-triggers mtime  -n all_demography
"""


rule all_demography:
    input: 
      sfs_prev = expand( "../results/demography/preview/prev_{spec}_on_{ref}.txt" , ref = "mirang", spec = "mirang" ),
      sfs_dir = expand( "../results/demography/sfs/{spec}_on_{ref}" , ref = "mirang", spec = "mirang" )

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