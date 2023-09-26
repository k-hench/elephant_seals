library(prismatic)
mm_to_inch <- 25.4
fwidth <- 178 / mm_to_inch
fhalfwidth <- 87 / mm_to_inch

# fnt_sel <- "Josefin sans"#"CMU Sans Serif"
fnt_sel <- "Arial"#"CMU Sans Serif"
fnt_sz <- 12 / ggplot2::.pt

specs <-  c("mirang", "mirleo") 
spec_names <- c(mirang = "M. angustirostris", mirleo = "M. leonina")
clrs <- c("#282828", "#db4a32") |> 
  purrr::set_names(specs)

clr_default <- c("#6495ED", "#777777")

clrs_n <- \(n){scales::colour_ramp(colors = clrs)((0:(n-1))/(n-1))}

clr_pheno <- RColorBrewer::brewer.pal(3, "Set1") |> 
  set_names(nm = c("worms", "control", "mirleo"))

clr_f_traj <- colorRampPalette(colors = c("black",
                                          clr_default[[2]], 
                                          clr_default[[1]], 
                                          clr_darken(clr_default[[1]])))

plt_lwd <- 0.15
clr_lines <- "black"

theme_ms <- \(..., lwd = .5){
  list(theme_bw(base_family = fnt_sel,
                ...),
       theme(panel.border = element_blank(),
             panel.grid = element_blank(),
             axis.line = element_line(colour = clr_lines,
                                      linewidth = lwd),
             axis.ticks = element_line(colour = clr_lines,
                                       linewidth = lwd),
             plot.subtitle = element_text(hjust = .5)))
}