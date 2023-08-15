library(tidyverse)

compl <- c(a = "T", A = "T", c = "G", C = "G", g = "C", G = "C", t = "A", T = "A")

data_anc <- read_tsv("results/ancestral_allele/mirang_anc52_snps.tsv.gz") |> 
  mutate(POS  = refPosition + 1)
data_snps <- read_tsv("results/ancestral_allele/snps_vcf.tsv.gz")

data_joined <- left_join(data_snps, data_anc,
                         by = c(`#CHROM` = "refSequence", POS = "POS")) |> 
  filter(complete.cases(refPosition)#,   # cases where mirang and anc are identical don't need updating
         #complete.cases(REF)
         )           # invariant cases in mirang don't need updating

chroms <- unique(data_snps$`#CHROM`)[1:50]

data_joined |> 
  filter(`#CHROM` == chroms[[2]]) |> 
  mutate(mir_compl = compl[mirang],
         anc_compl = compl[Anc52],
         mir_fw_org_check = mir_compl == REF,
         anc_fw_org_check = anc_compl == ALT,
         mir_fw_anti_check = 1 - (mir_compl == REF),
         anc_fw_anti_check = 1 - (anc_compl == ALT),
         mir_rev_org_check = mir_compl == ALT,
         anc_rev_org_check = anc_compl == REF,
         mir_rev_anti_check = 1 - (mir_compl == ALT),
         anc_rev_anti_check = 1 - (anc_compl == REF)) |>
  filter(mir_rev_org_check)
  summarise(across(ends_with("check"), sum, .names = "{.col}-count"),
            across(ends_with("check"), mean, .names = "{.col}-frq")) |> 
  pivot_longer(everything()) |> 
  separate(name, into = c("class", "type"), sep = "-") |> 
  separate(class, into = c("player", "dir", "logic", "tag"), sep = "_") |> 
  pivot_wider(names_from = c(type, logic), values_from = value) |> 
  mutate(frq_test = frq_org+frq_anti)
  
data_alt_separated <- data_joined |> 
    # filter(`#CHROM` == chroms[[2]]) |> 
    mutate(across(c(REF, ALT, mirang, Anc52), str_to_upper)) |> 
    separate(ALT, into = c("ALT_1", "ALT_2", "ALT_3", "ALT_4", "ALT_5"), sep = ",")

data_alt_separated |> 
    pivot_longer(starts_with("ALT"), values_to = "ALT") |> 
    select(-name) |> 
    filter(ALT %in% c("A", "C", "G", "T")) |> 
    mutate(# anc IS alt, swap alleles in vcf (default expected case)
           set_matches_c = ( REF == mirang & ALT == Anc52) | (REF == compl[mirang] & ALT == compl[Anc52]),
           # 
           alt_differs_c = !set_matches_c  & (( ALT != Anc52) | ( ALT != compl[Anc52])),
           ref_differs_c =  !set_matches_c  & (( REF != mirang) | (REF != compl[mirang]))#,
           #######
           # simple_matches_c = ( ALT == Anc52) | (ALT == compl[Anc52]),
           # # 
           # simple_alt_differs_c = !simple_matches_c  & (( ALT != Anc52) | ( ALT != compl[Anc52])),
           # simple_ref_differs_c =  !simple_matches_c  & (( REF != mirang) | (REF != compl[mirang]))
           ) |> 
  # filter(alt_differs_c)
    summarise(across(ends_with("_c"), sum, .names = "{.col}-pass"),
              across(ends_with("_c"), \(x){sum(1-x)}, .names = "{.col}-fail")) |> 
  pivot_longer(everything()) |> 
  arrange(name) |> 
  separate(name, into = c("check", "pass"), sep = "-") #|> 
  # group_by(check) |> 
  # summarise(sum(value))

data_alt_separated |> 
  pivot_longer(starts_with("ALT"), values_to = "ALT") |> 
  select(-name) |> 
  filter(ALT %in% c("A", "C", "G", "T")) |> 
  mutate(snp = str_c(`#CHROM`,"_", POS)) |> 
  pluck("snp") |> duplicated() |> sum()


data_alt_separated |> 
  filter( !( (REF == mirang) | (REF == compl[mirang])) ) |>
  # filter( !((REF == compl[mirang])) ) |>
  # filter( !((REF == mirang)) ) |>
  # filter(!(( REF != mirang) | (REF != compl[mirang])))|>
  # pivot_longer(starts_with("ALT"), values_to = "ALT") |> 
  # select(-name) |> 
  # filter(ALT %in% c("A", "C", "G", "T")) |> 
  # filter(!(( ALT != Anc52) | (ALT != compl[Anc52])))|> 
  mutate(comp = compl[mirang], 
         check = REF == comp)

# check if org or complementary is needed, based on REF > adjust Anc accordingly > THEN check if Anc == ALT or Anc %in% ALT...
data_check_anc <- data_joined |> 
  # filter(`#CHROM` == chroms[[2]])
  mutate(across(c(REF, ALT, mirang, Anc52), str_to_upper)) |> 
  mutate(ref_correct = case_when(
    REF == mirang ~ TRUE,
    REF == compl[mirang] ~ FALSE,
    .default = NA
  ),
  anc_corrected = if_else(ref_correct, Anc52, compl[Anc52]),
  anc_in_alt = map2_lgl(anc_corrected, ALT, grepl))# |> 
#  pluck("ref_correct") |> is.na() |> sum() 

data_check_anc |> filter(`#CHROM` == "NC_072356.1", POS == 44568)

# check how often REF is neither matching nor complement
sum(is.na(data_check_anc$ref_correct))
# check how often Anc in not listed in Alt
n_alt_non_match <- sum(1 - data_check_anc$anc_in_alt)

# fraction of diverging ancestral alleles not within the mirang population (relative to all diverging alleles)
n_alt_non_match / nrow(data_check_anc)
# fraction of diverging ancestral alleles not within the mirang population (relative to all variable sites)
n_alt_non_match / nrow(data_snps)
# fraction of variable sites where REF != Anc 
nrow(data_check_anc) / nrow(data_snps)

data_new_ref <- data_snps |> 
  left_join(data_check_anc |> 
              # remove diverging alleles that are not in the population
              filter(anc_in_alt) |>  
              select(`#CHROM`, POS, new_ref = anc_corrected),
            by = c("#CHROM", "POS")) |> 
  # if allele is not diverging in Anc, keep REF
  mutate(new_ref = if_else(is.na(new_ref), REF, new_ref))

data_new_ref |>
  filter(REF != new_ref) 
