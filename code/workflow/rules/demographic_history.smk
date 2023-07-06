"""
snakemake --rerun-triggers mtime  -n ../results/fastsimcoal/
"""


rule all_demography:
    input:  expand( "../results/demography/preview/prev_{spec}_on_{ref}.txt" , ref = "mirang", spec = ["mirang", "mirleo"] )

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