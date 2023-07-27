library(tidyverse)
library(here)
library(patchwork)
library(ggnewscale)

get_psl <- \(file){
  vroom::vroom(here::here("results","psl", file),
               delim = "\t",
               col_names = c("matches", "misMatches", "repMatches", "nCount",
                             "qNumInsert", "qBaseInsert", "tNumInsert",
                             "tBaseInsert", "strand", "qName", "qSize", "qStart",
                             "qEnd", "tName", "tSize", "tStart", "tEnd",
                             "blockCount")) |>
    select( tName, tStart, tEnd, qName, qStart, qEnd, strand ) |>
    mutate(tSize = abs(tEnd - tStart),
           qSize = abs(qEnd - qStart))
}

read_size <- \(genome = "", y_base = 0, skip = 0, order_by = "size", manual_order_18 = NULL, n_first = 18){
  read_tsv(here::here("results", "genomes", str_c(genome,".size")),
           col_names = c("chr", "size")) |>
    mutate(org_pos = row_number()) |>
    arrange(desc(size)) |>
    mutate(pre_manual = rank(-size),
           size_idx = if(order_by == "size"){
             row_number()
           } else {c(manual_order_18[[genome]], (n_first+1):length(chr))}) |>
    arrange(size_idx) |>
    mutate(end_with_skip = cumsum(size + skip),
           start = lag(end_with_skip, default =  0),
           end = start + size,
           mid = (start + end) /2,
           eo = row_number() %% 2,
           y_base = y_base,
           genome = genome,
           in_top = size_idx <= n_first)
}

data_check <- data_algn <- get_psl("slim_mirang_on_zalcal.psl.gz") 

data_check <- data_algn |> 
  filter(qName == "NC_072356.1")

data_check |> 
  group_by(tName, qName) |> 
  summarise(sum_length_t = sum(tSize),
            sum_length_q = sum(qSize),
            tSize = tSize[1],
            qSize = qSize[1],
            n_al = length(qSize)) |> 
  ggplot(aes(x = tName, y = sum_length_q/qSize)) +
  geom_point(aes(color = tName == "NC_045612.1"))


data_algn_x <- data_algn |> 
  filter(tName == "NC_045612.1")

data_summary <- data_algn |> 
  group_by(tName, qName) |> 
  summarise(sum_length_t = sum(tSize),
            sum_length_q = sum(qSize),
            n_al = length(qSize)) |> 
  left_join(read_size("mirang") |> select(qName = chr, qSize = size)) |> 
  left_join(read_size("zalcal") |> select(tName = chr, tSize = size)) |> 
  mutate(q_fract = sum_length_q / qSize,
         t_fract = sum_length_t / tSize) |> 
  ungroup()

data_summary_x <- data_summary |> 
  filter(qName %in% data_algn_x$qName) |> 
  group_by(qName) |> 
  mutate(q_lead = q_fract / sum(q_fract)) |> 
  ungroup() |> 
  group_by(qName) |> 
  arrange(qName) |> 
  mutate(q_x = (\(x,y){x[y == "NC_045612.1"]})(q_lead, tName),#q_lead[tName == "NC_045612.1"],
         q_max = max(q_lead,na.rm = TRUE),
         q_max_alt = (\(x,y){
           x_check <- max(x[y != "NC_045612.1"], na.rm = TRUE)
           if(x_check == -Inf){0}else{x_check}
           })(q_lead, tName),
         x_log_ratio = log10(q_x / q_max_alt),
         x_log_ratio = if_else(x_log_ratio == Inf, 3.5, x_log_ratio)) |> 
  ungroup()

data_summary_x |> 
  filter(qName == "NC_072356.1")

# p1 <- data_summary_x |> 
#   ggplot(aes(x = q_x, y = q_max_alt, color = log10(q_x / q_max_alt))) +
#   scico::scale_color_scico(limits = c(-3.5, 3.5), palette = "vanimo") +
#   geom_point(size = .2) +
#   theme_minimal()

treshold <- log10(2/1)

p1 <- ggplot(mapping = aes(x = q_x,
             y = q_max_alt)) +
  scico::scale_color_scico(limits = c(-3.5, 3.5), palette = "vanimo") +
  geom_point(data = data_summary_x |>
               filter(x_log_ratio > treshold),
             size = .2,
             aes(color = log10(q_x / q_max_alt))) +
  new_scale_color() +
  geom_point(data = data_summary_x |>
               filter(x_log_ratio < treshold),
             size = .2,
             aes(color = log10(q_x / q_max_alt))) + 
  scale_color_distiller(palette = "Greys", limits = c(-3.5, treshold)) +
  coord_equal(xlim = 0:1, ylim = 0:1) +
  theme_minimal() +
  theme(legend.position = "bottom")

p2 <- data_summary_x  |> 
  filter(x_log_ratio > treshold) |>
  arrange(-x_log_ratio) |> 
  mutate(q_f = fct_reorder(qName, -x_log_ratio),
         t_f = factor(tName)) |> 
  ggplot(aes(x = as.numeric(q_f),
             y = tName)) +
  geom_tile(aes(fill = x_log_ratio,
                color = tName == "NC_045612.1")) +
  scico::scale_fill_scico(limits = c(-3.5, 3.5), palette = "vanimo") +
  scale_color_manual(values = c(`TRUE` = "black", `FALSE` = "white")) +
  # coord_equal(xlim = c(0, 350)) +
  theme(legend.position = "bottom")

dir.create(here("results/genomes/sex_chrom"), showWarnings = FALSE)
data_summary_x  |> 
  filter(x_log_ratio > treshold) |> 
  filter(!duplicated(qName)) |> 
  mutate(start = 0) |> 
  select(chr = qName, start, end = qSize) |> 
  write_tsv(here("results/genomes/sex_chrom/mirang_sex_chrom.bed"))

read_size("mirang") |> 
   filter(chr == "NW_025578803.1")
