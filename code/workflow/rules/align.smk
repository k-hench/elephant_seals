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

ruleorder: unpack_genome > split_fasta
ruleorder: lastdb_index > gsalign_index
ruleorder: align_single_last > align_single_minimap > align_single_gsalign > align_split
ruleorder: lastdb_index > gsalign_index
ruleorder: align_split > align_single_last > align_single_minimap > align_single_gsalign 

rule align_all:
    message:
      """
      Aligning the elephant seal reference genomes to
      identify the sacffolds located on the X chromosome
      based on the hits on the zalcal scaffold `NC_045612.1`.
      """
    input:
      # the wildcards in psls will also handle definition of wildcards for fastas
      psls=expand("../results/psl/{species}.psl", species = ALIGN_SPECIES)

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
      temp('../results/genome/{refname}lastdb_index.done')
    params:
      indexBase='../data/genomes/{refname}',
      refSizeFile='../results/genome/{refname}.size',
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

rule gsalign_index:
    input:
      fasta='../data/genomes/{refname}.fa'
    output:
      # Create a dummy file that tells us indexing the ref has finished.
      # file is actually created in the shell: directive with "touch"
      # This rule will also create the ref's 2bit file, which may be
      # used later in net_to_axt (but isn't used at the time of writing)
      temp('../results/genome/{refname}lastdb_index.done')
    params:
      indexBase='../data/genomes/{refname}',
      refSizeFile='../results/genome/{refname}.size',
    log:
      'logs/{refname}_lastdbIndex_log.txt'
    threads: 2
    benchmark:
      'benchmark/{refname}-lastdb.txt'
    conda:
      'msa_align'
    shell:
      """
      faSize -detailed {input.fasta} > {params.refSizeFile} 2>{log} && GSAlign index {input.fasta} {params.indexBase} 2>{log} && touch {output} 2>{log}
      """

rule build_index:
    input:
      str(rules.lastdb_index.output).format(refname = ALIGN_REF),
      fastaFile="../data/genomes/{species}.fa",
      fastaRef='../data/genomes/{refname}.fa'.format(refname = ALIGN_REF)
    output:
      "../results/genome/{species}.size"
    params:
      indexBase='../data/genomes/{refname}'.format(refname = ALIGN_REF),
      speciesSizeFile='../results/genome/{species}.size',
      refNibDir='../results/genome/nib',
      refFastaFile='../data/genomes/{refname}.fa'.format(refname = ALIGN_REF),
      refNib2Bit='../results/genome/nib/{refname}.2bit'.format(refname = ALIGN_REF),
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

rule split_fasta:
    input:
      str(rules.lastdb_index.output).format(refname = ALIGN_REF),
      speciesSizeFile='../results/genome/{species}.size',
      fastaFile='../data/genomes/{species}.fa' 
    output:
      flat=temp('../data/genomes/{species}.fa.flat'),
      gdx=temp('../data/genomes/{species}.fa.gdx'),
      splitDummy=temp('../data/genomes/{species}.split')
    params:
      splitFastaN=1
    log:
      'logs/{species}_fasplit_log.txt'
    benchmark:
      'benchmark/{species}-split_bm.txt'
    conda:
      'msa_split'
    threads: 1
    shell:
      """
      pyfasta split -n {params.splitFastaN} {input.fastaFile} &>{log} && \
      touch {output.splitDummy}
      """

rule align_single_last:
    input:
      str(rules.lastdb_index.output).format(refname = ALIGN_REF),
      fastaFile="../data/genomes/{species}.fa",
      speciesSizeFile='../results/genome/{species}.size',
    output:
      name='../results/psl/{species}.psl'
    params:
      indexBase='../data/genomes/{refname}'.format(refname = ALIGN_REF),
      lastParams=LAST_PARAMS,
      lastSplitParams='',
    log:
      'logs/{species}_align_log.txt'
    benchmark:
      'benchmark/{species}-align_bm.txt'
    conda:
      'msa_align'
    threads: 1
    shell:
      # This shell will kick off for each fasta in the {genomedir}/fastas folder.  Each instance
      # of this rule gets 1 thread, but multiple lastal commands may be run, depending on the number of species
      # and the number of threads given on the command line.
      # NOTE: more threads means more memory used, and you could run out, so have to
      # temper the number of threads.
      # the file size from faSize is needed is the chain/net steps later as is the
      # ref .2bit file
      """
      lastal {params.lastParams} {params.indexBase} {input.fastaFile} {params.lastSplitParams} | maf-convert psl /dev/stdin 2>{log} 1>{output.name}
      """

rule align_single_minimap:
    input:
      str(rules.lastdb_index.output).format(refname = ALIGN_REF),
      fastaFile="../data/genomes/{species}.fa",
      speciesSizeFile='../results/genome/{species}.size',
    output:
      name='../results/psl/{species}.psl'
    params:
      minimap2Params = "-a -cx asm20",
      refFastaFile = '../data/genomes/{refname}.fa'.format(refname = ALIGN_REF),
    log:
      'logs/{species}_align_log.txt'
    benchmark:
      'benchmark/{species}-align_bm.txt'
    conda:
      'msa_align'
    threads: 1
    shell:
      # This shell will kick off for each fasta in the {genomedir}/fastas folder.  Each instance
      # of this rule gets 1 thread, but multiple lastal commands may be run, depending on the number of species
      # and the number of threads given on the command line.
      # NOTE: more threads means more memory used, and you could run out, so have to
      # temper the number of threads.
      # the file size from faSize is needed is the chain/net steps later as is the
      # ref .2bit file
      """
      minimap2 {params.minimap2Params} {params.refFastaFile} {input.fastaFile} 2>>{log} | samtools sort 2>>{log} | bamToPsl /dev/stdin {output.name} &>>{log}
      """

rule align_single_gsalign:
    input:
      str(rules.lastdb_index.output).format(refname = ALIGN_REF),
      fastaFile="../data/genomes/{species}.fa",
      speciesSizeFile='../results/genome/{species}.size',
    output:
      name='../results/psl/{species}.psl'
    params:
      indexBase='../data/genomes/{refname}'.format(refname = ALIGN_REF),
      speciesPath='../results/genome/{species}',
      gsalignParams="-sen -no_vcf",
    log:
      'logs/{species}_align_log.txt'
    benchmark:
      'benchmark/{species}-align_bm.txt'
    conda:
      'msa_align'
    threads: 1
    shell:
      # This shell will kick off for each fasta in the {genomedir}/fastas folder.  Each instance
      # of this rule gets 1 thread, but multiple lastal commands may be run, depending on the number of species
      # and the number of threads given on the command line.
      # NOTE: more threads means more memory used, and you could run out, so have to
      # temper the number of threads.
      # the file size from faSize is needed is the chain/net steps later as is the
      # ref .2bit file
      """
      GSAlign {params.gsalignParams} -r {params.refFastaFile} -q {input.fastaFile} -o {params.speciesPath} -i {params.indexBase} &>>{log} && \
      sed -i 's/^s qry\./s /' {params.speciesPath}.maf && sed -i 's/^s ref\./s /' {params.speciesPath}.maf && \
      maf-convert psl {params.speciesPath}.maf 2>>{log} 1>{output.name}
      """

rule align_split:
    input:
      splitFa= "../data/genomes/{species}.fa",
      splitDummy='../data/genomes/{species}.split',
      speciesSizeFile='../results/genome/{species}.size'
    output:
      psl='../results/psl/{species}.psl',
      cmd=temp('../results/genome/{species}.cmd'),
      cmdcomp=temp('../results/genome/{species}.cmd.completed'),
      splitMaf=temp('../results/genome/{species}.maf')
    params:
      indexBase='../data/genomes/{refname}'.format(refname = ALIGN_REF),
      refName=ALIGN_REF,
      splitDir='../results/genome/',
      speciesPath='../results/genome/{species}',
      lastParams=LAST_PARAMS,
      lastSplitParams='',
      splitFa='../data/genomes/{species}.fa'
    log:
      'logs/{species}_lastAlign_split_log.txt'
    benchmark:
      'benchmark/{species}_lastAlign_split_bm.txt'
    conda:
      'msa_align'
    threads: 1
    shell:
      # This script will align split fasta files to the reference using parafly parallelization
      # It uses some potentially unsafe globbing and rm
      # These should be replaced with expand() inputs by eliminating 0 padding from split fasta names
      """
      lastal {params.lastParams} {params.indexBase} data/genomes/${input.splitFa} > {output.splitMaf} &> {log}
      sed '30,${{/^#/d;}}' {output.splitMaf} | maf-sort /dev/stdin {params.lastSplitParams} | maf-convert psl /dev/stdin | awk '$9!="++"' > {output.psl}
      """
