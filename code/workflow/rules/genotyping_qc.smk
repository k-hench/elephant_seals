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

// **Allelic imbalance**
//
// More straight forward alternative for contamination-check through allelic imbalance (based on [workshop by Joana Meier and Mark Ravinet](https://speciationgenomics.github.io/allelicBalance/))
//
// The idea here is to use the coverage ratio of the different alleles at heterozygous sites to screen for PCR duplications and cross-sample contamination. 
// To be able to differentiate between those two potential souces of allelic imbalance, we compare the whole SNP set whit the subset of SNPs that have a minor allele count (MAC) of 1 within the data set.
//
// So we start by creating the MAC filtered subset for the genotypes.
//
// ```groovy
process filter_mac1 {
  label 'mac1_c_pop'
  time { 1.h * task.attempt }
  memory '5. GB'

  input:
  tuple val( spec ), file( filt ), file( biallele )

  output:
  tuple file( "${biallele[0]}" ), file( "*_mac1.vcf.gz" )

  script:
  """
  FL_BASE=\$(echo ${biallele[0]} | sed 's/.vcf.gz//')
  vcftools \
    --gzvcf ${biallele[0]} \
    --max-mac 1 \
    --mac 1 \
    --recode \
    --stdout |  \
    bgzip > \${FL_BASE}_mac1.vcf.gz
  """
}
// ```
//
// For each individual in the `.vcf` file, we then print the coverage at both alleles for every heterozygous SNP.
//
// ```groovy
process export_het_ind {
  label 'exp_het_c_pop'
  time { 5.h * task.attempt }
  memory '25. GB'
  scratch "${baseDir}/../../tmp"

  input:
  tuple path( vcf )

  output:
  tuple file( vcf ), file( "*_hets_*.csv" ), file( "*.tmp.indv" )

  script:
  """
  r=\$(echo ${vcf} | sed 's/.vcf.gz//')
  # Creates individual file
  zgrep "#CHROM" ${vcf} | \
    cut -f 10- | \
    sed 's/\\t/\\n/g' > \$r".tmp.indv"
  
  # Writes the AD fields from all heterozygous genotypes with minor allele count i into a file hets.i
  n=\$(wc -l \$r".tmp.indv" | tr -s " " | cut -f 1 -d " ")
  c=1
  for i in \$(cat \$r".tmp.indv")
  do
    echo -ne "\$i | \$c of \$n"\\\\r
          vcftools \
            --gzvcf ${vcf} \
            --mac 1 \
            --max-mac 1 \
            --indv \$i \
            --recode \
            --stdout |\
          grep -v "^#" | \
          grep "0[/|]1" |\
          awk '{split(\$9,a,":");for(i=1;i<=10;i++){if(match(a[i],"AD")){adidx=i}};for(i=10;i<=NF;i++){if(match(\$i,"0[/|]1")==1){split(\$i,b,":");print b[adidx]}}}' \
          > \$r"_hets_"\$i".csv"
    ((c=c+1))
  done
  """
}
// ```
//
// The allele coverage is then binned in 2 dimensions (Allele-count for "A", and allele count for "B") for every sample.
//
// ```groovy
process bin_hets {
  publishDir "${params.outdir}/qc/allelic_imbalance", mode: 'copy', pattern: ".tsv"
  label "bin_hets_c_renv_v1"
  time { 3.h * task.attempt }
  memory '20. GB'

  input:
  tuple file( vcf ), file( hets ), file( inds )

  output:
  tuple file( vcf ), file( "*hetIndStats_freq2d.tsv" ), file( "*hetIndStats_d.tsv" )

  script:
  """
  #!/usr/bin/env Rscript
  library(tidyverse)
  base <- str_remove("${vcf}", ".vcf.gz");
  s <- base; # used to be random number
  xmaxx <- 0;

  freq_2d_all <- tibble(x = c(), y = c(), het_GT_count= c())
  d_all <- tibble(V1= c(), V2= c(), max= c(), min= c(), minreadprop= c(), depth= c(), ind = c())
  inds <- as.character(read.table("${inds}")\$V1);

  for(i in inds){
    if(file.exists(paste(base,"_hets_",i,".csv",sep=""))){if(file.info(paste(base,"_hets_",i,".csv",sep=""))\$size>0){
      d<-read.table(paste(base,"_hets_",i,".csv",sep=""),sep=",");
      d<-cbind(d,max=apply(d,1,max),min=apply(d,1,min));
      d<-cbind(d,minreadprop=d\$min/(d\$min+d\$max),depth=d\$min+d\$max);
      x.bin<-seq(floor(min(d\$depth,na.rm=T)),ceiling(max(d\$depth,na.rm=T)),length=ceiling(max(d\$depth,na.rm=T))-floor(min(d\$depth,na.rm=T))+1);
      y.bin<-seq(floor(min(d\$min,na.rm=T)),ceiling(max(d\$min,na.rm=T)),length=ceiling(max(d\$min,na.rm=T))-floor(min(d\$min,na.rm=T))+1);
      freq<-as.data.frame(table(findInterval(d\$depth,x.bin),findInterval(d\$min,y.bin)));
      freq[,1] <- as.numeric(as.character(freq[,1]));
      freq[,2] <- as.numeric(as.character(freq[,2]));
      freq2D<-matrix(0,nrow=length(x.bin),ncol=length(y.bin));
      freq2D[cbind(freq[,1], freq[,2])] <- freq[,3];

      freq_i <- freq2D |> 
        as_tibble() |> 
        rename_with(.fn = \\(x){str_remove(x, "V")}) |> 
        mutate(x = row_number()) |> 
        pivot_longer(cols = -x, names_to = "y",names_transform = as.numeric,values_to = 'het_GT_count') |> 
        mutate(ind = i)
      d_i <- d |> 
        as_tibble() |> 
        mutate(ind = i)
      
      freq_2d_all <- freq_2d_all |> 
        bind_rows(freq_i)
      
      d_all <- d_all |> 
        bind_rows(d_i)
    }}
  };

  freq_2d_all |> write_tsv(paste(base,"hetIndStats_freq2d.tsv",sep="."))
  d_all |> write_tsv(paste(base,"hetIndStats_d.tsv",sep="."))
  """
}
// ```
//
// Finally a 2D heatmap of the allele counts, a histogram and the density of the allelic ratios is plotted.
//
// ```groovy
process plot_allelic_imbalance {
  publishDir "${params.outdir}/fig/qc/allelic_imbalance", mode: 'copy'
  label "bin_hets_c_renv_v1"
  time { 1.h * task.attempt }
  memory '12. GB'

  input:
  tuple file( vcf ), file( ind_freq ), file( ind_d )

  output:
  file( "*.pdf" )

  script:
  """
  #!/usr/bin/env Rscript
  library(tidyverse)

   base <- str_remove("${vcf}", ".vcf.gz");

  library(tidyverse)
  library(prismatic)
  library(patchwork)
  library(glue)

  d <- read_tsv(glue("{base}.hetIndStats_d.tsv"))
  freq2D <- read_tsv(glue("{base}.hetIndStats_freq2d.tsv"))

  clrs <- c("#fcde9c", "#faa476", "#f0746e",
            "#e34f6f", "#dc3977", "#b9257a", "#7c1d6f")
  clr_pick <- clrs[[2]]
  clr_high <- 'black'
  
  seq_depth <- 30
  p1 <- freq2D |> 
    ggplot(aes(x =x ,y = y, fill = log10(het_GT_count))) +
    geom_raster() +
    geom_abline(data = tibble(a = c(.5, .3,.2)),
                aes(slope = a, intercept = 0, linetype = factor(a)),
                alpha = .5, linewidth = .4) +
    facet_wrap(ind~., nrow = 2)+
    # scale_fill_viridis_c(option = "A", na.value = 'transparent') +
    scale_fill_gradientn(colours = clrs, na.value = 'transparent',
                        breaks = (c(0:6)/2), labels = \\(x){sprintf('%.1f',10^x)}) +
    scale_linetype_manual(values = c(`0.5` = 1, `0.3` = 2, `0.2` = 3),
                          guide = "none")+
    scale_x_continuous(limits = c(-1, 3.5 * seq_depth)) +
    scale_y_continuous(limits = c(-1, 1.75 * seq_depth)) +
    guides(fill = guide_colorsteps(title.position = 'top',
                                  title = "log10(hetrozygous genotype count)",
                                  barwidth = unit(.75, "npc"),
                                  barheight = unit(5,"pt"))) +
      labs(x = "DP", y = "Reads minor alleles") +
    coord_cartesian(xlim = c(0, 3 * seq_depth),
                    ylim = c(0, 1.5 * seq_depth))

  p2 <- d |> 
    ggplot(aes(x = minreadprop)) +
    geom_histogram(binwidth = 0.02, boundary = 0,
                  color = clr_pick,
                  fill = clr_lighten(clr_pick))+
    geom_vline(data = tibble(a = c(.5, .3,.2)),
                aes(xintercept =  a, linetype = factor(a)),
              alpha = .5, linewidth = .4) +
    scale_linetype_manual(values = c(`0.5` = 1, `0.3` = 2, `0.2` = 3),
                          guide = "none") +
    facet_wrap(ind~., nrow = 2)+
    labs(x = "Proportion minor allele reads") +
    coord_cartesian(xlim = c(0,0.5))

  p3 <- d |> 
    ggplot(aes(x = minreadprop)) +
    #geom_rect(data = tibble(x = Inf, ind = 105407),
    #          inherit.aes = FALSE,
    #          aes(xmin =-x, xmax = x, ymin =-x, ymax = x),
    #          color = clr_high, fill = clr_alpha(clr_high,.05)) +
    geom_density(adjust = 0.4,
                  color = clr_pick,
                  fill = clr_lighten(clr_pick))+
    geom_vline(data = tibble(a = c(.5, .3,.2)),
                aes(xintercept =  a, linetype = factor(a)),
              alpha = .5, linewidth = .4) +
    scale_linetype_manual(values = c(`0.5` = 1, `0.3` = 2, `0.2` = 3),
                          guide = "none") +
    facet_wrap(ind~., nrow = 2)+
    labs(x = "Proportion minor allele reads") +
    coord_cartesian(xlim = c(0,0.5))


  p_out <- p1 / p2 / p3 + 
    plot_annotation(subtitle = base) +
    plot_layout(guides = 'collect') &
    theme_minimal() &
    theme(legend.position = "bottom",
          plot.subtitle = element_text(hjust = .5))

  ggsave(plot = p_out,
        filename = glue("{base}_het_stats_all_ind.pdf"),
        width = 16,
        height = 12,
        device = cairo_pdf)
  """
}
'''
