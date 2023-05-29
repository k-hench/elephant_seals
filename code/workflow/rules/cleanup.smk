UNPACKED_FA=[ 'mirang', 'mirleo', 'zalcal' ]

rule cleanup_all:
    input: 
      file = expand("../data/genomes/{spec}.fa",  spec = UNPACKED_FA)
    shell:
      """
      rm {input.file}
      """