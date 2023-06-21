'''
snakemake -n -R  align_all
snakemake --dag -R  align_all | dot -Tsvg > ../results/img/control/dag_align.svg
'''

from snakemake import available_cpu_count
from psutil import virtual_memory

available_mem_gb = lambda: '%dG' % (virtual_memory().available >> 30)
ALIGN_REF = "zalcal"
ALIGN_SPECIES = ['mirang', 'mirleo']

# split fasta index numbers
splitN = 1
#padN = len(str(splitN))
#Nlist = list(range(0, splitN))
#padList = [str(item).zfill(padN) for item in Nlist]

# Prevent evaluation of lastParams as None
LAST_PARAMS = "-m 10 -j 3 -u 1 -p HOXD70"

minimap2Params: "-a -cx asm20"
gsalignParams: "-sen -no_vcf"


rule align_all:
    message:
      """
      Aligning the elephant seal reference genomes to
      identify the sacffolds located on the X chromosome
      based on the hits on the zalcal scaffold `NC_045612.1`.
      """
    input:
      # the wildcards in psls will also handle definition of wildcards for fastas
      psls=expand("../results/psl/slim_{species}_on_{ref}.psl.gz", species = ALIGN_SPECIES, ref = ALIGN_REF )

rule unpack_genome:
    input:
      gz ='../data/genomes/{refname}.fa.gz'
    output:
      fasta = "../data/genomes/{refname}.fa"
    conda:
      'msa_align'
    shell:
      """
      zcat {input.gz} > {output.fasta}
      """

rule lastdb_index:
    input:
      fasta='../data/genomes/{refname}.fa'
    output:
      # Create a dummy file that tells us indexing the ref has finished.
      # file is actually created in the shell: directive with "touch"
      # This rule will also create the ref's 2bit file, which may be
      # used later in net_to_axt (but isn't used at the time of writing)
      temp('../results/genomes/{refname}lastdb_index.done')
    params:
      indexBase='../data/genomes/{refname}',
      refSizeFile='../results/genomes/{refname}.size',
    log:
      'logs/align/{refname}_lastdbIndex_log.txt'
    threads: 2
    benchmark:
      'benchmark/{refname}-lastdb.txt'
    conda:
      'msa_align'
    container: "docker://khench/msa_envs:v0.1"
    shell:
      """
      faSize -detailed {input.fasta} > {params.refSizeFile} 2>{log} && lastdb -R 10 -u YASS -c {params.indexBase} {input.fasta} 2>{log} && touch {output} 2>{log}
      """

rule build_index:
    input:
      str(rules.lastdb_index.output).format(refname = ALIGN_REF),
      fastaFile="../data/genomes/{species}.fa",
      fastaRef='../data/genomes/{refname}.fa'.format(refname = ALIGN_REF)
    output:
      "../results/genomes/{species}.size"
    params:
      indexBase='../data/genomes/{refname}'.format(refname = ALIGN_REF),
      speciesSizeFile='../results/genomes/{species}.size',
      refNibDir='../results/genomes/nib',
      refFastaFile='../data/genomes/{refname}.fa'.format(refname = ALIGN_REF),
      refNib2Bit='../results/genomes/nib/{refname}.2bit'.format(refname = ALIGN_REF),
    log:
      'logs/{species}_index_log.txt'
    benchmark:
      'benchmark/{species}-index_bm.txt'
    conda:
      'msa_ucsc'
    threads: 1
    shell:
      """
      mkdir -p {params.refNibDir} && \
      faToTwoBit {params.refFastaFile} {params.refNib2Bit} && \
      faSize -detailed {input.fastaFile} > {output}
      """

rule align_single_last:
    input:
      str(rules.lastdb_index.output).format( refname = ALIGN_REF ),
      fastaFile = "../data/genomes/{species}.fa.gz",
      speciesSizeFile = '../results/genomes/{species}.size',
    output:
      maf = '../results/maf/{species}_on_{ref}.maf.gz'
    params:
      indexBase = '../data/genomes/{ref}'.format( ref = ALIGN_REF ),
      lastParams = LAST_PARAMS,
      mafBase = '../results/maf/{species}_on_{ref}.maf'
    log:
      'logs/align/{species}_on_{ref}_align.log'
    conda:
      'msa_align'
    threads: 1
    shell:
      """
      lastal {params.lastParams} {params.indexBase} {input.fastaFile} 2>{log} 1>{params.mafBase}
      gzip {params.mafBase}
      """

rule maf_to_psl:
    input:
      maf = '../results/maf/{species}_on_{ref}.maf.gz'
    output:
      psl = '../results/psl/{species}_on_{ref}.psl.gz'
    params:
      pslBase = '../results/psl/{species}_on_{ref}.psl'
    log:
      'logs/psl/{species}_on_{ref}_psl.log'
    conda:
      'msa_align'
    threads: 1
    shell:
      """
      zcat {input.maf} | maf-convert psl 2>{log} 1>{params.pslBase}
      gzip {params.pslBase}
      """

rule slim_psl:
    input: '../results/psl/{species}_on_{ref}.psl.gz'
    output: '../results/psl/slim_{species}_on_{ref}.psl.gz'
    shell:
      """
      zcat {input} | cut -f 1-18  | gzip > {output}
      """