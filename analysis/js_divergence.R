#### Setup ####

library(philentropy)
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

hap_freq <- read_csv("data/processed/poolseq/haplotype_freq_200kb.csv") %>%
  # remove samples that were resequenced (they are highly correlated)
  filter(!(sample %in% c("HN-directional-B", "HN-stabilising-A", "HN-stabilising-B", "LN-directional-A", "LN-random-A", "LN-stabilising-B"))) %>%
  mutate(sample = str_replace(sample, "_reseq", ""))

# major allele frequencies
maf <- hap_freq %>%
  group_by(sample, chrom, pos) %>%
  filter(freq == max(freq)) %>%
  ungroup() %>%
  group_by(sample, nitrate, selection, rep, chrom, pos, freq) %>%
  summarise(acc = paste(unique(acc), collapse = " ")) %>%
  ungroup() %>%
  mutate(nacc = str_count(acc, " ") + 1)



# Divergence --------------------------------------------------------------

jsdist <- hap_freq %>%
  select(sample, chrom, pos, acc, freq) %>%
  pivot_wider(names_from = "sample", values_from = "freq") %>%
  group_nest(chrom, pos) %>%
  mutate(jsd = map(data, function(i){
    i %>%
      column_to_rownames("acc") %>%
      as.matrix() %>%
      t() %>%
      distance(method = "jensen-shannon", unit = "log2",
               use.row.names = TRUE, as.dist.obj = TRUE) %>%
      broom::tidy()
  })) %>%
  unnest(jsd) %>%
  select(-data)

# filter to retain only those compared with respective control
jsdist %>%
  # same nitrate treatment
  filter(
    (str_detect(item1, "HN") & str_detect(item2, "HN")) |
      str_detect(item1, "LN") & str_detect(item2, "LN")
    ) %>%
  # same replicate
  filter(
    (str_detect(item1, "A") & str_detect(item2, "A")) |
      (str_detect(item1, "B") & str_detect(item2, "B")) |
      (str_detect(item1, "C") & str_detect(item2, "C"))
  ) %>%
  # compared to random control
  filter(
    (str_detect(item1, "directional") & str_detect(item2, "random")) |
      (str_detect(item2, "directional") & str_detect(item1, "random"))
  ) %>%
  ggplot(aes(pos, distance)) +
  geom_line(aes(group = interaction(item1, item2))) +
  facet_grid(item2 ~ chrom, scales = "free_x", space = "free_x")


# compare MajAF between selected and controls
maf %>%
  group_by(nitrate, rep, chrom, pos) %>%
  summarise(diff = freq[selection == "directional"] - freq[selection == "random"]) %>%
  ungroup() %>%
  ggplot(aes(pos, diff)) +
  geom_line(size = 1, colour = "steelblue") +
  geom_hline(yintercept = 0, lty = 2) +
  facet_grid(nitrate + rep ~ chrom, scales = "free_x", space = "free_x")



# Allele frequency --------------------------------------------------------

# allele frequency distributions
hap_freq %>%
  ggplot(aes(freq)) +
  stat_ecdf(aes(group = sample, colour = selection)) +
  facet_wrap(~ acc) +
  geom_vline(xintercept = 1/19, colour = "brown")

# major allele frequency distributions
maf %>%
  ggplot(aes(freq)) +
  stat_ecdf(aes(group = sample, colour = selection)) +
  facet_wrap( ~ nitrate)

# replicate frequency spectrum
maf %>%
  filter(freq > 0.3) %>%
  group_by(chrom, pos, nitrate, selection) %>%
  summarise(nsample = n_distinct(sample)) %>%
  ungroup() %>%
  count(nitrate, selection, nsample) %>%
  ggplot(aes(nsample, n)) +
  geom_line(aes(colour = selection), size = 1) +
  geom_point(pch = 21, fill = "white", size = 3) +
  facet_wrap(~ nitrate) +
  scale_x_continuous(breaks = 1:3)

# consistent alelle
hap_freq %>%
  filter(freq > 0.3) %>%
  distinct(sample, chrom, pos, nitrate, selection, acc) %>%
  group_by(chrom, pos, nitrate, selection) %>%
  summarise(nsample = n_distinct(sample), nacc = n_distinct(acc)) %>%
  ungroup() %>%
  count(nitrate, selection, nsample, nacc) %>%
  ggplot(aes(nsample, n, colour = selection)) +
  geom_line(size = 1) +
  geom_point(pch = 21, fill = "white", size = 3) +
  facet_grid(nacc ~ nitrate) +
  scale_x_continuous(breaks = 1:3)

hap_freq %>%
  filter(freq > 0.3) %>%
  distinct(sample, chrom, pos, nitrate, selection) %>%
  group_by(chrom, pos, nitrate, selection) %>%
  summarise(nsample = n_distinct(sample)) %>%
  ungroup() %>%
  count(nitrate, selection, nsample) %>%
  ggplot(aes(nsample, n, colour = selection)) +
  geom_line(size = 1) +
  geom_point(pch = 21, fill = "white", size = 3) +
  facet_grid( ~ nitrate) +
  scale_x_continuous(breaks = 1:3)

hap_freq %>%
  filter(selection == "random") %>%
  group_by(sample, acc) %>%
  summarise(q99 = quantile(freq, 0.99)) %>%
  group_by(acc) %>%
  summarise(q99 = median(q99)) %>%
  full_join(hap_freq, by = "acc") %>%
  filter(freq > q99) %>%
  distinct(sample, chrom, pos, nitrate, selection) %>%
  group_by(chrom, pos, nitrate, selection) %>%
  summarise(nsample = n_distinct(sample)) %>%
  ungroup() %>%
  count(nitrate, selection, nsample) %>%
  ggplot(aes(nsample, n, colour = selection)) +
  geom_line(size = 1) +
  geom_point(pch = 21, fill = "white", size = 3) +
  facet_grid( ~ nitrate) +
  scale_x_continuous(breaks = 1:3)


# trying to set a MajAF threshold based on the selected populations
thr <- maf %>%
  filter(selection == "directional") %>%
  group_by(sample) %>%
  summarise(q99 = quantile(freq, 0.99)) %>%
  ungroup() %>%
  summarise(q99 = median(q99)) %>%
  pull(q99)

maf %>%
  mutate(sig = ifelse(freq > thr, freq, NA)) %>%
  filter(selection != "stabilising") %>%
  ggplot(aes(pos, freq)) +
  geom_line(aes(colour = selection)) +
  geom_point(aes(y = sig)) +
  facet_grid(nitrate + rep ~ chrom, scales = "free_x", space = "free_x")


# setting a threshold on individual populations
maf %>%
  group_by(sample) %>%
  mutate(sig = ifelse(freq > quantile(freq, 0.99), freq, NA)) %>%
  ungroup() %>%
  filter(selection != "stabilising") %>%
  ggplot(aes(pos, freq)) +
  geom_line(aes(colour = selection)) +
  geom_point(aes(y = sig, colour = selection)) +
  facet_grid(nitrate + rep ~ chrom, scales = "free_x", space = "free_x")



# Filter windows ----------------------------------------------------------

temp <- hap_freq %>%
  select(sample, chrom, pos, acc, freq) %>%
  group_by(chrom, pos) %>%
  filter(any(freq > 0.2)) %>%
  ungroup() %>%
  pivot_wider(names_from = "sample", values_from = "freq") %>%
  group_nest(chrom, pos) %>%
  mutate(jsd = map(data, function(i){
    i %>%
      column_to_rownames("acc") %>%
      as.matrix() %>%
      t() %>%
      distance(method = "jensen-shannon", unit = "log2",
               use.row.names = TRUE, as.dist.obj = TRUE) %>%
      broom::tidy()
  })) %>%
  unnest(jsd)

temp <- temp %>%
  filter(!str_detect(item1, "stabilising") & !str_detect(item2, "stabilising")) %>%
  filter(!str_detect(item1, "reseq") & !str_detect(item2, "reseq")) %>%
  filter((str_detect(item1, "HN") & str_detect(item2, "HN")) |
           str_detect(item1, "LN") & str_detect(item2, "LN")) %>%
  mutate(comparison = case_when(
    str_detect(item1, "directional") &
      str_detect(item2, "directional") ~
      "directional vs directonal",
    (str_detect(item1, "directional") &
       str_detect(item2, "random")) | (str_detect(item2, "directional") &
                                         str_detect(item1, "random")) ~
      "directional vs random",
    str_detect(item1, "random") &
      str_detect(item2, "random") ~
      "random vs random",
    TRUE ~  NA_character_))


temp %>%
  ggplot(aes(distance)) +
  geom_density(aes(colour = comparison, group = interaction(item1, item2)))

