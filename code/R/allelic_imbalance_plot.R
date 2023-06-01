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