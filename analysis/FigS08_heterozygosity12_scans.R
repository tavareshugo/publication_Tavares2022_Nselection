#
# Fig S07
#

#### Setup ####

library(tidyverse)
library(patchwork)

# Change ggplot2 defaults
theme_set(theme_bw() +
            theme(panel.grid = element_blank(),
                  text = element_text(size = 10),
                  legend.background = element_blank()))

scale_colour_continuous <- scale_colour_viridis_c
scale_colour_discrete <- function(palette = "Dark2", ...) scale_colour_brewer(palette = palette, ...)
scale_fill_continuous <- scale_fill_viridis_c
scale_fill_discrete <- function(palette = "Dark2", ...) scale_fill_brewer(palette = palette, ...)


#### Read data ####

source("./analysis/data_processing/read_poolseq.R")

# Define sweep that will be focused on in Chr1 region
target_sweep1 <- candidate_sweeps %>%
  filter(sample == "LN-directional-C" & chrom == 1 & statistic == "H12")

# Define sweep that will be focused on in Chr5 region
target_sweep5 <- candidate_sweeps %>%
  filter(nitrate == "HN" & chrom == 5 & statistic == "H12")




#### panel A ####

p1 <- candidate_sweeps %>%
  filter(selection != "stabilising") %>%
  mutate(statistic = str_remove(statistic, "H"),
         selection = factor(selection, levels = c("random", "directional", "stabilising"))) %>%
  ggplot() +
  geom_segment(aes(x = start/1e6, xend = end/1e6, y = statistic, yend = statistic,
                   colour = selection),
               size = 2) +
  geom_point(data = pool200, aes(pos/1e6, 1 - hap_hom12), alpha = 0) +
  geom_point(data = bind_rows(target_sweep1, target_sweep5),
             aes(x = pos/1e6, y = 0.3), shape = 17, size = 3) +
  geom_rect(data = centromeres,
            aes(xmin = start/1e6 - 0.5, xmax = end/1e6 + 0.5, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.5) +
  facet_grid(nitrate + rep ~ chrom, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = seq(0, 30, 10)) +
  labs(x = "Mb",
       y = expression(paste(italic(H[e]), " statistic"))) +
  theme_bw() +
  theme(legend.position = "none")


#### panel B ####

# get frequencies around the H12 peak on LN - B
sweep_freqs <- target_sweep1 %>%
  select(chrom, pos, selection) %>%
  left_join(pool200, by = c("chrom", "pos", "selection")) %>%
  filter(nitrate == "LN") %>%
  select(sample, nitrate, selection, rep, hap_hom1, freqs_csv) %>%
  drop_na(freqs_csv) %>%
  mutate(freqs = map(freqs_csv, read_csv, col_names = c("accession", "frequency"))) %>%
  select(-freqs_csv) %>%
  unnest(freqs)

# check top 2 accessions in each population
sweep_freqs %>%
  group_by(sample) %>%
  top_n(2, frequency)

# Haplotype frequency spectrum for each candidate sweep
p2 <- sweep_freqs %>%
  ggplot(aes(frequency, str_to_sentence(accession))) +
  geom_col() +
  facet_grid( ~ paste(nitrate, rep, sep = " - ")) +
  labs(x = "Frequency", y = "Accession",
       title = "Chr 1 at 9.4Mb")


#### panel C ####

# get frequencies around the H12 peak on LN - B
sweep_freqs <- target_sweep5 %>%
  select(chrom, pos, selection) %>%
  left_join(pool200, by = c("chrom", "pos", "selection")) %>%
  filter(nitrate == "HN") %>%
  select(sample, nitrate, selection, rep, hap_hom1, freqs_csv) %>%
  drop_na(freqs_csv) %>%
  mutate(freqs = map(freqs_csv, read_csv, col_names = c("accession", "frequency"))) %>%
  select(-freqs_csv) %>%
  unnest(freqs)

# check top 2 accessions in each population
sweep_freqs %>%
  group_by(sample) %>%
  top_n(2, frequency) %>%
  arrange(sample)

# Haplotype frequency spectrum for each candidate sweep
p3 <- sweep_freqs %>%
  ggplot(aes(frequency, str_to_sentence(accession))) +
  geom_col() +
  facet_grid( ~ paste(nitrate, rep, sep = " - ")) +
  labs(x = "Frequency", y = "Accession", 
       title = "Chr 5 at 8Mb")

# save graph
p1 + p2 + p3 +
  plot_layout(ncol = 1, heights = c(2.5, 1, 1)) +
  plot_annotation(tag_levels = "A")
ggsave("./figures/FigS08.pdf", width = 7.5, height = 10)

# visualise the whole frequency spectrum across chromosome
pool200 %>%
  filter(chrom == 5 & nitrate == "HN" & selection == "directional") %>%
  drop_na(freqs_csv) %>%
  mutate(freqs = map(freqs_csv, read_csv, col_names = c("accession", "frequency"))) %>%
  select(-freqs_csv) %>%
  unnest() %>%
  filter(frequency > 0.05) %>%
  arrange(frequency) %>%
  ggplot() +
  geom_point(aes(x = pos/1e6, y = accession, colour = frequency), size = 2) +
  geom_rect(data = centromeres %>% filter(chrom == 1),
            aes(xmin = start/1e6 - 0.5, xmax = end/1e6 + 0.5, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.5) +
  facet_grid(rows = vars(rep)) +
  scale_colour_gradient2(low = "white", mid = "steelblue", high = "brown", midpoint = 0.4) +
  geom_vline(xintercept = 8000000/1e6)
