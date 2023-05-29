"""
snakemake -c 1 --use-conda geno_prep
snakemake --dag -R  geno_prep | dot -Tsvg > ../results/img/control/dag_geno_prep.svg
"""

rule geno_prep:
    input: '../results/checkpoints/done_round_0.txt'
