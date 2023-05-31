"""
"""

rule all_gt_qc:
    message:
      """
      Check sequencing quality (`fastqc`), mapping coverage
      (`bamcov` and `bamtools`),
      inter-specific contamination (`fastq-screen`),
      allelic imbalance for intra-specific contamination.
      """


'''
process run_fastqc {
  publishDir "${params.outdir}/qc/", mode: 'copy'
  tag "${sample}.${info.lane}"
  label "Q_fastqc_c_qc"
  memory '15. GB'
  time '2.5h'

  input:
  tuple val( spec ), val( sample ), val( fw ), val( rv ),  val( info )

  output:
  file( "fastqc/*" )

  script:
  """
  mkdir -p fastqc
  fastqc ${params.seq_dir}/${fw} ${params.seq_dir}/${rv} -o fastqc/
  """
}

process subsample_fastq {
  label "Q_def_subsample_c_qc"
  memory '8. GB'
  tag "${sample}.${info.lane}"

  input:
  tuple val( spec ), val( sample ), val( fw ), val( rv ),  val( info )

  output:
  tuple val( "${sample}" ), file( "*.1.fq.gz" ), file( "*.2.fq.gz" )

  script:
  """
  seqtk sample -s 42 ${params.seq_dir}/${fw} 100000 | bgzip > ${sample}.${info.lane}.1.fq.gz
  seqtk sample -s 42 ${params.seq_dir}/${rv} 100000 | bgzip > ${sample}.${info.lane}.2.fq.gz
  """
}

process run_fastq_screen {
  publishDir "${params.outdir}/qc/fastq_screen", mode: 'copy'
  tag "${sample}"
  label "fast_screen_c_qc"
  time { 20.m * task.attempt }
  memory '3. GB'
  clusterOptions '-V -cwd -P fair_share -l idle=1 -l si_flag=1 -pe multislot 7'

  input:
  tuple val( sample ), file( fw ), file( rv )

  output:
  tuple file( "*_screen.html" ), file( "*_screen.txt" )

  script:
  """
  sed "s=<fq_screen_dir>=${params.base}/data/fq_screen_db/=g" \
    ${params.base}/data/fq_screen_db/fastq_screen.conf > fastq_screen.conf

  fastq_screen --conf fastq_screen.conf --threads 7 --aligner bowtie2  ${fw} ${rv}
  """
}

process run_bamcov {
  publishDir "${params.outdir}/qc/coverage", mode: 'copy'
  label 'Q_def_bamcov_c_qc'
  memory '40. GB'
  tag "${sample}.${info.lane}"

  input:
  tuple val( sample ), val( spec ), val( lane ), val( info ), file( bam ), file( bai ), file( tsv )

  output:
  file( "*_coverage.tsv.gz" )

  script:
  """
  bamcov ${bam} -o ${sample}.${info.lane}_on_${params.reference}_coverage.tsv
  gzip ${sample}.${info.lane}_on_${params.reference}_coverage.tsv
  """
}

process run_bamcov_mask {
  publishDir "${params.outdir}/qc/coverage/masks", mode: 'copy'
  label 'Q_def_bamcov_c_pop'
  memory '40. GB'
  tag "${sample}.${info.lane}"

  input:
  tuple val( sample ), val( spec ), val( lane ), val( info ), file( bam ), file( bai ), file( tsv )

  output:
  file( "*_covmask.bed.gz" )

  script:
  """
  samtools view -b ${bam} | \
    genomeCoverageBed -ibam stdin -bg | \
    gzip > ${sample}.${info.lane}_on_${params.reference}_covmask.bed.gz
  """
}

process run_bamtools {
  publishDir "${params.outdir}/qc/bamstats", mode: 'copy'
  label 'Q_def_bamtools_c_qc'
  memory '20. GB'
  tag "${sample}.${info.lane}"

  input:
  tuple val( sample ), val( spec ), val( lane ), val( info ), file( bam ), file( bai ), file( tsv )

  output:
  file( "*_on_${params.reference}.bamstats" )

  script:
  """
  bamtools stats -in ${bam} | \
    grep -v "*" > ${sample}_${info.lane}_on_${params.reference}.bamstats
  """
}
'''