dat = [ "a", "b" ]

rule all_examples:
  input: expand("{bd}/results/data_{d}.tsv", d = dat, bd = base_dir)

rule copy_data:
  input: "{bd}/data/data_{d}.tsv"
  output: "{bd}/results/data_{d}.tsv"
  log: "{bd}/code/logs/copy_{d}.log"
  shell:
    """
    echo {wildcards.d} > {log}
    cp {input} {output}
    """