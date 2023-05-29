library(tidyverse)
library(ape)
library(ggtree)
library(patchwork)
library(here)

data <- read_tsv(here("data", "carnivora_ref_genomes.tsv")) |>
  mutate(`Organism Name` = `Organism Name` |> 
           str_replace("Canis lupus familiaris", "Canis lupus") |> 
           str_replace("Taxidea taxus jeffersonii", "Taxidea taxus") |> 
           str_replace("Enhydra lutris kenyoni", "Enhydra lutris") |> 
           str_replace("Mustela putorius furo", "Mustela putorius") |> 
           str_replace("Odobenus rosmarus divergens", "Odobenus rosmarus") |> 
           str_replace("Ursus thibetanus thibetanus", "Ursus thibetanus"),
         repo = c(GCF = "refseq", GCA = "genbank")[str_sub(`Assembly Accession`,1,3)] |> 
           factor(levels = c("refseq", "genbank"))) |>
  arrange(`Organism Name`,
          as.numeric(repo), 
          -as.numeric(`Assembly Submission Date`))
#  filter(!duplicated(`Organism Name`, fromLast = TRUE))
tree <- read.tree(here("data", "carnivora_list.nwk"))

data_tree <- ggtree(tree, layout = "fan",
                    open.angle = 180)$data |> 
  groupClade(.node = 119) |> 
  mutate(spec_group = if_else(label == "Mirounga_angustirostris", "focal",
                              if_else(label %in% c("Mirounga_leonina", "Leptonychotes_weddellii"),
                                      "outgroup", "other")),
         label = str_replace(label, "_", " "))

name_to_spec_lab <- \(nm){
  str_c(str_sub(nm, 1, 3),
        str_extract(nm, "_[a-z]{3}") |> str_remove("_")) |> 
  str_to_lower()
}

table_export <- data_tree |> 
  arrange(label) |> 
  filter(grepl( "[a-z]", label)) |>
  dplyr::select(`Organism Name` = label, spec_group) %>%
  left_join(data |>
              dplyr::select(`Organism Name`,
                            `Assembly Accession`,
                            `Assembly Submission Date`,
                            repo), . ) |>
  filter(!is.na(spec_group)) |> 
  # filter(!duplicated(`Organism Name`)) # 67
  group_by(`Organism Name`) |> 
  mutate(genome_alternative = row_number()) |> 
  ungroup() |> 
  mutate(`Organism Name` = str_replace_all(`Organism Name`, " ", "_") |> 
           str_to_lower(),
         genome_version = str_c(`Organism Name`, "_", genome_alternative),
         spec = name_to_spec_lab(`Organism Name`)) |>
  filter(genome_alternative == 1, `Organism Name` != "Mirounga_leonina") |>
  set_names(nm = function(str){str_to_lower(str) |> str_replace_all(" ", "_")}) |>
  dplyr::select(spec, everything())

p <- data_tree  |> 
  mutate(organism_name = str_replace_all(label, " ", "_") |> 
           str_to_lower() ) |> 
  left_join(table_export) |> 
  ggtree(layout = "fan",
         aes(color = spec_group),
         open.angle = 15, size = .4) +
  geom_tiplab(aes(x = x+ 2.5,
                  label = str_c("[ ", spec, " ]")),
                  family = "Ubuntu mono",
              size = 3) +
  geom_tiplab(aes(x = x+ 24.5,
                  label = str_replace(label, "_", " ")),
              size = 3) +
  scale_x_continuous(limits = c(0, 120)) +
  scale_color_manual(values = c(focal = "#4f7ca6", outgroup = "black", other = "gray50"),
                     guide = 'none')

ggsave(filename = here("results", "img", "reference_phylogeny.pdf"),
       plot = p,
       width = 9, height = 9, device = cairo_pdf)

# check if any tree labels are missing from the export table
# sum(!str_to_lower(tree$tip.label) %in% table_export$organism_name)

tree_short_labels <- tree
tree_short_labels$tip.label <- name_to_spec_lab(tree$tip.label)
tree_short_labels |> write.tree(here("data", "carnivora_short_labels.nwk"))

# (ggtree(tree) + geom_tiplab() + coord_cartesian(xlim = c(0, 67))) +
# (ggtree(tree_short_labels) + geom_tiplab(hjust = 1) + coord_cartesian(xlim = c(60, -5)))

table_export |> 
  write_tsv(here("data", "carnivora_genome_and_timetree.tsv"))

table_export |> 
  bind_rows(tibble(spec = "galgal",
                   organism_name = "gallus_gallus",
                   assembly_accession = "GCF_016699485.2",
                   assembly_submission_date = as_date("2021-01-19"),
                   repo = "refseq",
                   spec_group = "other",
                   genome_alternative = 1,
                   genome_version = "gallus_gallus_1")) |>
  write_tsv(here("data", "carnivora_and_chicken.tsv"))