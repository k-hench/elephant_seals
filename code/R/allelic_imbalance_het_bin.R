 #!/usr/bin/env Rscript
  library(tidyverse)
  base <- str_remove("${vcf}", ".vcf.gz");
  s <- base; # used to be random number
  xmaxx <- 0;

  freq_2d_all <- tibble(x = c(), y = c(), het_GT_count= c())
  d_all <- tibble(V1= c(), V2= c(), max= c(), min= c(), minreadprop= c(), depth= c(), ind = c())
  inds <- as.character(read.table("${inds}")$V1);

  for(i in inds){
    if(file.exists(paste(base,"_hets_",i,".csv",sep=""))){if(file.info(paste(base,"_hets_",i,".csv",sep=""))$size>0){
      d<-read.table(paste(base,"_hets_",i,".csv",sep=""),sep=",");
      d<-cbind(d,max=apply(d,1,max),min=apply(d,1,min));
      d<-cbind(d,minreadprop=d$min/(d$min+d$max),depth=d$min+d$max);
      x.bin<-seq(floor(min(d$depth,na.rm=T)),ceiling(max(d$depth,na.rm=T)),length=ceiling(max(d$depth,na.rm=T))-floor(min(d$depth,na.rm=T))+1);
      y.bin<-seq(floor(min(d$min,na.rm=T)),ceiling(max(d$min,na.rm=T)),length=ceiling(max(d$min,na.rm=T))-floor(min(d$min,na.rm=T))+1);
      freq<-as.data.frame(table(findInterval(d$depth,x.bin),findInterval(d$min,y.bin)));
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