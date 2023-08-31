library(ggplot2)
library(purrr)
library(prismatic)
# the default font family
fnt_sel <- "Arial"

# the default plotting color as well as the secondary plotting color
clr_default <- c("#6495ED", "#777777")

# the color scheme for the phenotypes
clr_pheno <- RColorBrewer::brewer.pal(8, "Set1")[c(1,3:5,7:8)] |>
  color() |> 
  set_names(nm = c("Worms", "Bacteria", "Malnutrition"," Trauma", "Protozoa", "Congenital defect"))

# the color scheme for the load types
clr_load <- c("#000000", "gray80", "gray35", "gray60") |> 
  color() |> 
  set_names(nm = "total", "masked", "fixed", "expressed")

clr_f_traj <- colorRampPalette(colors = c("gray90",
                                          clr_default[[1]], 
                                          clr_darken(clr_default[[1]],.85)))

clr_lines <- "black"

theme_ms <- \(..., lwd = .3){
  list(theme_bw(base_family = fnt_sel,
                ...),
       theme(panel.border = element_blank(),
             panel.grid = element_blank(),
             axis.line = element_line(colour = clr_lines,
                                      linewidth = lwd),
             axis.ticks = element_line(colour = clr_lines,
                                       linewidth = lwd),
             axis.text = element_text(color = "black"),
             plot.subtitle = element_text(hjust = .5),
             strip.background = element_blank()))
}
