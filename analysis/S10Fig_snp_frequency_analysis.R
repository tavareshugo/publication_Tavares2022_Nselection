#
# Fig S10
#

library(tidyverse)
library(valr)
library(patchwork)

theme_set(theme_bw() +
            theme(panel.grid = element_blank(),
                  text = element_text(size = 10),
                  legend.background = element_blank()))

scale_colour_continuous <- scale_colour_viridis_c
scale_colour_discrete <- function(palette = "Dark2", ...) scale_colour_brewer(palette = palette, ...)
scale_fill_continuous <- scale_fill_viridis_c
scale_fill_discrete <- function(palette = "Dark2", ...) scale_fill_brewer(palette = palette, ...)


# Read data ------

# read snape frequencies
snape200 <- read_csv("./data/processed/poolseq/snp_freq_200kb.csv") %>%
  mutate(selection = factor(selection, levels = c("random", "directional", "stabilising")))

# Read centromere region as annotated in TAIR9 release
# ftp://ftp.arabidopsis.org/home/tair/Genes/TAIR9_genome_release/TAIR9_gff3/Assembly_GFF
# see https://www.biostars.org/p/18782/
centromeres <- read_csv("./data/external/tair9_centromeres.csv",
                        col_types = cols())

# read heterozygosity
source("./analysis/data_processing/read_poolseq.R")

# annotate each window according to whether it falls within 1MB of the centromere
snape200 <- snape200 %>%
  distinct(chrom, start = win_start, end = win_end) %>%
  bed_intersect(centromeres %>% 
                mutate(centromere = TRUE, 
                       start = start - 1e6, 
                       end = end + 1e6)) %>%
  select(chrom, start = start.x, end = end.x, centromere = centromere.y) %>%
  full_join(snape200, by = c("chrom", "start" = "win_start", "end" = "win_end")) %>%
  mutate(centromere = replace_na(centromere, FALSE))


# Define outliying windows -----

# using 99% percentile
# centromeres are excluded
snape200 <- snape200 %>%
  group_by(selection) %>%
  mutate(outlier = ifelse(afc_mean > quantile(afc_mean, 0.99), afc_mean, NA)) %>%
  ungroup() %>%
  mutate(outlier = ifelse(centromere, NA, outlier))

# define peaks (for plotting)
peaks <- snape200 %>%
  filter(!is.na(outlier)) %>%
  group_by(nitrate, rep, selection, chrom) %>%
  # merge those within 400kb
  bed_cluster(max_dist = 400e3) %>%
  group_by(nitrate, rep, selection, chrom, .id) %>%
  summarise(start = min(start), 
            end = max(end), 
            win_mid = win_mid[which(afc_mean == min(afc_mean))],
            afc_mean = min(afc_mean)) %>%
  group_by(selection) %>%
  arrange(chrom, start) %>%
  mutate(peak_id = 1:n()) %>%
  ungroup() %>%
  select(-.id)


# Exploratory ------

# mean |AFC|
snape200 %>%
  filter(selection != "stabilising") %>%
  ggplot() +
  geom_line(aes(win_mid/1e6, afc_mean, colour = selection)) +
  geom_point(data = peaks, shape = 25, size = 1,
             aes(win_mid/1e6, afc_mean + 0.05, fill = selection, colour = selection)) +
  geom_rect(data = centromeres,
            aes(xmin = start/1e6 - 0.5, xmax = end/1e6 + 0.5, ymin = -Inf, ymax = Inf), fill = "grey", alpha = 0.5) +
  facet_grid(nitrate + rep ~ chrom, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = seq(0, 30, 10)) +
  labs(x = "Mb", colour = "selection: ") +
  theme(legend.position = "top") +
  coord_cartesian(ylim = c(0, 0.3))

# Proportion of |AFC| outliers per window
snape200 %>%
  filter(selection != "stabilising") %>%
  ggplot() +
  geom_line(aes(win_mid/1e6, afc_outlier_mean, colour = selection)) +
  geom_rect(data = centromeres,
            aes(xmin = start/1e6 - 0.5, xmax = end/1e6 + 0.5, ymin = -Inf, ymax = Inf), fill = "grey", alpha = 0.5) +
  facet_grid(nitrate + rep ~ chrom, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = seq(0, 30, 10)) +
  labs(x = "Mb", colour = "selection: ") +
  theme(legend.position = "top") +
  coord_cartesian(ylim = c(0, 0.3)) + 
  labs(y = "Fraction |AFC| outliers")

# these two are highly correlated anyway
snape200 %>%
  filter(selection != "stabilising" & !centromere) %>%
  ggplot(aes(afc_mean, afc_outlier_mean, colour = !is.na(outlier))) +
  geom_point()

# Correlation between He and |AFC|
# we gave this figure in a reviewer reply
pool200 %>%
  select(-freqs_csv) %>%
  inner_join(snape200,
             by = c("nitrate", "selection", "rep", "chrom", "pos" = "win_mid")) %>%
  filter(!centromere) %>%
  mutate(outlier = case_when(!is.na(outlier) & hap_het < threshold$hap_hom1_q1 ~ "Both",
                             !is.na(outlier) ~ "SNP-based", 
                             hap_het < threshold$hap_hom1_q1 ~ "He-based", 
                             TRUE ~ NA_character_)) |> 
  mutate(outlier = fct_infreq(outlier)) |> 
  arrange((as.numeric(outlier))) |> 
  ggplot(aes(hap_het, afc_mean)) +
  geom_point(aes(colour = outlier)) +
  facet_grid(selection ~ nitrate + rep) +
  labs(x = "Heterozygosity", y = "|AFC|", colour = "Outlier") +
  scale_colour_brewer(palette = "Dark2", na.value = "grey50")


snape200 %>% 
  filter(!is.na(outlier)) |> 
  arrange(chrom, win_mid)

# one window which is only visible in snape data
snape200 %>%
  filter(selection == "directional", nitrate == "LN", rep == "A") %>%
  filter(chrom == 2, win_mid > 10e6) %>%
  arrange(desc(afc_mean))

# let's check what it's doing in the heterozygosity
pool200 %>%
  filter(selection == "directional", nitrate == "LN", rep == "A", pos == 13100000, chrom == 2) %>%
  mutate(freqs = map(freqs_csv, ~ read.csv(text = .x, header = FALSE, col.names = c("acc", "freq")))) %>%
  pull(freqs) %>%
  bind_rows() %>%
  ggplot(aes(freq, acc)) +
  geom_col()

# why didn't we pick this up with He123?
pool200 %>%
  filter(selection != "stabilising") %>%
  ggplot() +
  geom_line(aes(pos/1e6, 1-hap_hom123, colour = selection)) +
  geom_hline(yintercept = threshold$hap_hom1_q1, linetype = 2, colour = "grey48") +
  geom_rect(data = centromeres,
            aes(xmin = start/1e6 - 0.5, xmax = end/1e6 + 0.5, ymin = -Inf, ymax = Inf), fill = "grey", alpha = 0.5) +
  geom_segment(data = candidate_sweeps,
               aes(x = start/1e6, xend = end/1e6, y = 0.3, yend = 0.3), size = 3) +
  facet_grid(nitrate + rep ~ chrom, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = seq(0, 30, 10)) +
  labs(x = "Mb", y = "Heterozygosity", colour = "selection: ") +
  theme(legend.position = "top")


# Paper Figure -----

# genome scan
p1 <- snape200 %>%
  filter(selection != "stabilising") %>%
  ggplot() +
  geom_line(aes(win_mid/1e6, afc_mean, colour = selection)) +
  geom_point(data = peaks, shape = 25, size = 2,
             aes(win_mid/1e6, afc_mean + 0.05, fill = selection, colour = selection), 
             show.legend = FALSE) +
  geom_rect(data = centromeres,
            aes(xmin = start/1e6 - 0.5, xmax = end/1e6 + 0.5, ymin = -Inf, ymax = Inf), fill = "grey", alpha = 0.5) +
  facet_grid(nitrate + rep ~ chrom, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = seq(0, 30, 10)) +
  labs(x = "Mb", colour = "selection: ") +
  theme(legend.position = "top") +
  coord_cartesian(ylim = c(0, 0.3)) +
  labs(x = "Mb", y = "|AFC|")
  
# allele frequency for two new peaks
chr2_peak <- pool200 %>%
  filter(selection == "directional", 
         nitrate == "LN", pos == 13100000, chrom == 2) %>%
  mutate(freqs = map(freqs_csv, ~ read.csv(text = .x, header = FALSE, col.names = c("acc", "freq")))) %>%
  unnest(freqs) %>%
  ggplot(aes(freq, acc)) +
  geom_col() +
  facet_grid(~ rep) +
  labs(x = "Frequency", y = "Accession", 
       title = "Chr 2 at 13MB")

chr3_peak <- pool200 %>%
  filter(selection == "directional", 
         nitrate == "LN", pos == 2200000, chrom == 3) %>%
  mutate(freqs = map(freqs_csv, ~ read.csv(text = .x, header = FALSE, col.names = c("acc", "freq")))) %>%
  unnest(freqs) %>%
  ggplot(aes(freq, acc)) +
  geom_col() +
  facet_grid(~ rep) +
  labs(x = "Frequency", y = "Accession", 
       title = "Chr 3 at 2.2MB")

# save figure
p1 + chr2_peak + chr3_peak +
  plot_layout(ncol = 1, heights = c(2.5, 1, 1)) +
  plot_annotation(tag_levels = "A")
ggsave("./figures/S10Fig.pdf", width = 7.5, height = 10)
