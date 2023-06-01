"""
snakemake -n -R all_gt_qc
snakemake --jobs 5 --use-singularity --singularity-args "--bind $CDATA" --use-conda -R all_gt_qc

snakemake --jobs 60 \
  --latency-wait 30 \
  -p \
  --default-resources mem_mb=51200 threads=1 \
  --use-singularity \
  --singularity-args "--bind $CDATA" \
  --use-conda \
  --cluster '
    qsub \
      -V -cwd \
      -P fair_share \
      -l idle=1 \
      -l si_flag=1 \
      -pe multislot {threads} \
      -l vf={resources.mem_mb}' \
  --jn job_gt.{name}.{jobid}.sh \
  -R all_gt_qc && mv job_gt.* logs/
"""

rule all_gt_qc:
    message:
      """
      Check sequencing quality (`fastqc`), mapping coverage
      (`bamcov` and `bamtools`),
      inter-specific contamination (`fastq-screen`),
      allelic imbalance for intra-specific contamination.
      """
    input:
      by_sample_ln = expand( ["../results/checkpoints/fastqc/{sample_ln}.check",
                              "../results/qc/fastq_screen/{sample_ln}_fw_screen.txt"],
                              sample_ln = SAMPLES_LN ),
      by_sample = expand(["../results/qc/coverage/{sample_id}_on_{ref}_coverage.tsv.gz",
                          "../results/qc/coverage/masks/{sample_id}_on_{ref}_covmask.bed.gz",
                          "../results/qc/bamstats/{sample_id}_on_{ref}.bamstats"],
                         sample_id = SAMPLES, ref = GATK_REF[0])
      # GATK_REF[0] <- subset to mirang for now for disc-usage

rule fastqc:
    input:
      fq_fw = lambda wc: "../data/raw_sequences/" + get_sample_info(wc, what = "file_fw"),
      fq_rv = lambda wc: "../data/raw_sequences/" + get_sample_info(wc, what = "file_rv")
    output: 
      check = touch( "../results/checkpoints/fastqc/{sample_ln}.check" )
    benchmark:
      'benchmark/qc/fastqc_{sample_ln}.tsv'
    resources:
      mem_mb=10240
    container: c_qc
    shell:
      """
      fastqc {input.fq_fw} {input.fq_rv} -o ../results/qc/fastqc/
      """

rule subsample_fastq:
    input:
      fq_fw = lambda wc: "../data/raw_sequences/" + get_sample_info(wc, what = "file_fw"),
      fq_rv = lambda wc: "../data/raw_sequences/" + get_sample_info(wc, what = "file_rv")
    output:
      fq_fw = temp( "../data/raw_sequences/subsets/{sample_ln}_fw.fq.gz" ),
      fq_rv = temp( "../data/raw_sequences/subsets/{sample_ln}_rv.fq.gz" )
    benchmark:
      'benchmark/qc/subsample_{sample_ln}.tsv'
    resources:
      mem_mb=8192
    container: c_qc
    shell:
      """
      seqtk sample -s 42 {input.fq_fw} 100000 | bgzip > {output.fq_fw}
      seqtk sample -s 42 {input.fq_rv} 100000 | bgzip > {output.fq_fw}
      """

rule fastq_screen:
    input:
      fq_fw = "../data/raw_sequences/subsets/{sample_ln}_fw.fq.gz",
      fq_rv = "../data/raw_sequences/subsets/{sample_ln}_rv.fq.gz"
    output:
      fw = "../results/qc/fastq_screen/{sample_ln}_fw_screen.txt",
      rv = "../results/qc/fastq_screen/{sample_ln}_rv_screen.txt"
    benchmark:
      'benchmark/qc/fq_screen_{sample_ln}.tsv'
    resources:
      mem_mb=8192
    threads: 7
    container: c_qc
    shell:
      """
      fastq_screen \
        --conf fastq_screen.conf \
        --threads 7 \
        --outdir ../results/qc/fastq_screen/ \ 
        --aligner bowtie2 {input.fq_fw} {input.fq_rv}
      """

rule merge_bam_by_sample:
    input:
      bams = lambda wc: "../results/mapped_bams/" + gather_sample_entries(wc, what = "sample_ln") + "_on_{ref}.dedup.bam"
    output:
      single_bam = '../results/mapped_bams/combined_per_sample/{sample_id}_on_{ref}.bam'
    benchmark:
      'benchmark/qc/merge_bam_{sample_id}_on_{ref}.tsv'
    resources:
      mem_mb=51200
    container: c_gatk
    shell:
      """
      BAMS=$( echo {input.bams} | sed "s/\[//g; s/\]//g; s/,//g; s/'//g" )
      samtools merge  -o {output.single_bam} $BAMS
      gatk --java-options "-Xmx45g" BuildBamIndex -I {output.single_bam}
      """

rule bamcov:
    input:
      bam = "../results/mapped_bams/combined_per_sample/{sample_id}_on_{ref}.bam"
    output:
      cov = "../results/qc/coverage/{sample_id}_on_{ref}_coverage.tsv.gz"
    benchmark:
      'benchmark/qc/bamcov_{sample_id}_on_{ref}.tsv'
    params:
      unzipped = "../results/qc/coverage/{sample_id}_on_{ref}_coverage.tsv"
    resources:
      mem_mb=51200
    container: c_qc
    shell:
      """
      bamcov {input.bam} -o {params.unzipped}
      gzip {params.unzipped}
      """

rule bamcov_mask:
    input:
      bam = "../results/mapped_bams/combined_per_sample/{sample_id}_on_{ref}.bam"
    output:
      bed = "../results/qc/coverage/masks/{sample_id}_on_{ref}_covmask.bed.gz"
    benchmark:
      'benchmark/qc/covmask_{sample_id}_on_{ref}.tsv'
    resources:
      mem_mb=40960
    container: c_popgen
    shell:
      """
      samtools view -b {input.bam} | \
        genomeCoverageBed -ibam stdin -bg | \
        gzip > {output.bed}
      """

rule bamtools:
    input:
      bam = "../results/mapped_bams/combined_per_sample/{sample_id}_on_{ref}.bam"
    output:
      stats = "../results/qc/bamstats/{sample_id}_on_{ref}.bamstats"
    benchmark:
      'benchmark/qc/bamtools_{sample_id}_on_{ref}.tsv'
    resources:
      mem_mb=40960
    container: c_qc
    script:
      """
      bamtools stats -in {input.bam} | \
        grep -v "*" > {output.stats}
      """

'''
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
