library(ggplot2)
library(purrr)
library(prismatic)

# the default font family
fnt_sel <- "Arial"

# the default plotting color as well as the secondary plotting color
clr_default <- c("#253741", "#7198AD")

# the color scheme for the phenotypes
clr_pheno <- c("#E41A1C", "#4DAF4A", "#984EA3", "#FF7F00", "#DEC102", "#E4A1A0") |> 
  color() |> 
  set_names(nm = c("Worms", "Bacteria", "Malnutrition", "Trauma", "Protozoa", "Congenital defect"))

# the color scheme for the load types
clr_load <- c("#000000", "#ADD2E8", "#253741", "#7198AD") |> 
  color() |> 
  set_names(nm = "total", "masked", "fixed", "expressed")

clr_f_traj <- colorRampPalette(colors = c("gray90",
                                          clr_default[[1]], 
                                          clr_darken(clr_default[[1]],.85)))

point_sz <- 2

clr_lines <- "black"

theme_ms <- \(..., lwd = .2){
  list(theme_bw(base_family = fnt_sel,
                ...),
       theme(panel.border = element_blank(),
             panel.grid = element_blank(),
             axis.line = element_line(colour = clr_lines,
                                      linewidth = lwd),
             axis.ticks = element_line(colour = clr_lines,
                                       linewidth = lwd),
             axis.text = element_text(color = "black"),
             axis.text.x = element_text(color = "black"),
             axis.text.y = element_text(color = "black"),
             plot.subtitle = element_text(hjust = .5),
             strip.background = element_blank()))
}

spec_names <- c(mirang = "M. angustirostris", mirleo = "M. leonina")
