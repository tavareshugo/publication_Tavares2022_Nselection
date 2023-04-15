#
# Fig 05
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


#### Panel A ####

# Define chromosome sizes
chrom_sizes <- pool200 %>%
  group_by(chrom) %>%
  summarise(min = min(pos/1e6), max = max(pos/1e6)) %>%
  ungroup()

# QTL location from de Jong 2019 Table S1
qtl <- tibble(chrom = factor(c(1, 5, 1, 3, 2, 5)),
              pos = c(25735094, 26691233, 26274750, 5910414, 18811694, 5452597),
              nitrate = c("HN", "HN", "LN", "LN", "Plas", "Plas"))

# This figure has to be edited in inkscape
p1 <- candidate_sweeps %>%
  filter(statistic %in% c("H1", "H12") & selection == "directional") %>%
  mutate(peak_id = ifelse(statistic == "H12", NA, peak_id),
         pop_id = paste0(toupper(nitrate), " - ", toupper(rep)),
         colour_id = paste0(toupper(nitrate), " - ", statistic)) %>%
  ggplot() +
  geom_linerange(data = chrom_sizes, aes(x = as.character(chrom), ymin = min, ymax = max),
                 size = 1) +
  geom_linerange(aes(x = pop_id,
                     ymin = start/1e6, ymax = end/1e6,
                     colour = colour_id),
                 size = 3, position = position_dodge(width = 0.7)) +
  geom_text(aes(x = pop_id, y = (start + end)/2/1e6, label = peak_id),
            vjust = 2.2, size = 2.5) +
  geom_linerange(data = centromeres,
                 aes(x = as.character(chrom), ymin = start/1e6-0.5, ymax = end/1e6+0.5),
                 size = 3, colour = "grey", alpha = 0.8) +
  facet_grid(chrom ~ ., scales = "free", space = "free") +
  geom_label(data = mutate(qtl, chrom = as.character(chrom)), 
             aes(x = chrom, y = pos/1e6, label = nitrate),
             label.padding = unit(0.1, "lines"), size = 2.5) +
  scale_y_continuous(breaks = seq(0, 30, 5)) +
  coord_flip() +
  scale_colour_manual(values = c("LN - H1" = "#084594",
                                 "LN - H12" = "#bdd7e7",
                                 "HN - H1" = "#b2182b",
                                 "HN - H12" = "#fcae91")) +
  theme(strip.background = element_blank(), legend.position = "none",
        panel.grid.major.y = element_line(colour = "lightgrey")) +
  labs(x = "", y = "Mb", tag = "A", colour = "Nitrate - Statistic")


#### Panel B ####

# Haplotype frequency spectrum for each candidate sweep
sweep_freqs <- candidate_sweeps %>%
  filter(selection == "directional" & statistic == "H1") %>%
  distinct(sample, chrom, pos, peak_id, selection) %>%
  left_join(pool200, by = c("chrom", "pos", "selection")) %>%
  mutate(selected = sample.x == sample.y) %>%
  select(sample = sample.y, peak_id, nitrate, selection, rep, hap_hom1, freqs_csv, selected) %>%
  drop_na(freqs_csv) %>%
  mutate(freqs = map(freqs_csv, read_csv, col_names = c("accession", "frequency"))) %>%
  select(-freqs_csv) %>%
  unnest(cols = c(freqs))

# HN - A peaks
p2.1 <- sweep_freqs %>%
  filter(peak_id %in% c(4, 9) & sample == "HN-directional-A") %>%
  mutate(accession = str_to_sentence(accession)) %>%
  ggplot(aes(accession, frequency)) +
  geom_col() +
  facet_grid(cols = vars(peak_id)) +
  labs(x = "Accession", y = "Frequency", title = "HN - A", tag = "B") +
  scale_y_continuous(limits = c(0, 0.85)) +
  coord_flip()

# LN - A peaks
p2.2 <- sweep_freqs %>%
  filter(peak_id %in% c(3, 8) & sample == "LN-directional-A_reseq") %>%
  mutate(accession = str_to_sentence(accession)) %>%
  ggplot(aes(accession, frequency)) +
  geom_col() +
  facet_grid(cols = vars(peak_id)) +
  labs(x = "Accession", y = "Frequency", title = "LN - A") +
  scale_y_continuous(limits = c(0, 0.85)) +
  coord_flip()

# LN - B peaks
p2.3 <- sweep_freqs %>%
  filter(peak_id %in% c(1, 12) & sample == "LN-directional-B") %>%
  mutate(accession = str_to_sentence(accession)) %>%
  ggplot(aes(accession, frequency)) +
  geom_col() +
  facet_grid(cols = vars(peak_id)) +
  labs(x = "Accession", y = "Frequency", title = "LN - B") +
  scale_y_continuous(limits = c(0, 0.85)) +
  coord_flip()

# LN - C peaks
p2.4 <- sweep_freqs %>%
  filter(peak_id %in% c(7, 10) & sample == "LN-directional-C") %>%
  mutate(accession = str_to_sentence(accession)) %>%
  ggplot(aes(accession, frequency)) +
  geom_col() +
  facet_grid(cols = vars(peak_id)) +
  labs(x = "Accession", y = "Frequency", title = "LN - C") +
  scale_y_continuous(limits = c(0, 0.85)) +
  coord_flip()



p1 /
{p2.1 + p2.2 + p2.3 + p2.4 + plot_layout(ncol = 2, nrow = 2)} +
  plot_layout(ncol = 1, heights = c(1,1.5))
ggsave("./figures/Fig05.pdf", width = 7.5, height = 10)
