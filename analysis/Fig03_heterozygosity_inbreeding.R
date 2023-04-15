#
# Fig 03
#

#### Setup ####

library(tidyverse)
library(patchwork)
library(here)

# ensure working directory is project root
setwd(here())

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

# read phenotype data
phen_sum <- readRDS("./data/processed/phenotypes/phenotypes_summarised.rds")
# source("./analysis/data_processing/read_phenotypes.R")

# heterozygosity scans at 200Kb windows
pool200 <- read_csv("./data/processed/poolseq/haplotype_diversity_200kb.csv")
pool200 <- pool200 %>%
  mutate(selection = factor(selection, levels = c("random", "directional", "stabilising")),
         hap_het = 1 - hap_hom1) %>%
  select(sample, nitrate, selection, rep, chrom, pos, hap_het)

# read starting population scans at 200kb intervals
startpop <- read_csv("./data/external/founder_genotypes/startpop_window200kb.csv")
startpop <- startpop %>%
  mutate(sample = "founder population") %>%
  select(sample, chrom, pos, hap_het)
# startpop200 <- read_csv("./data_processed/startpop/startpop_window200kb.csv",
#                         col_types = "iiiiidididdddddddddddddddddd")

# put these two together
pool200 <- bind_rows(pool200, startpop)


#### make graphs ####

# boxplot of heterozygosity
h0 <- 1 - (19 * ((1/19)^2)) # starting heterozygosity
h10 <- h0 * (1 - 1/(2*40))^10 # theoretical gen10 heterozygosity
p1 <- pool200 %>%
  group_by(sample) %>%
  mutate(med = median(hap_het)) %>%
  ungroup() %>%
  mutate(sample = reorder(sample, med)) %>%
  ggplot(aes(sample, hap_het, fill = selection, group = sample)) +
  geom_hline(yintercept = c(h0, h10),
             linetype = "dashed", colour = "grey48") +
  geom_boxplot(position = position_dodge2(preserve = "total")) +
  geom_point(stat = "summary", fun = "mean") +
  labs(x = "Sample", y = "Heterozygosity", tag = "A") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.7))


# inbreeding coef across generations
p2 <- phen_sum %>%
  ggplot(aes(generation, f_mean)) +
  geom_line(aes(colour = selection)) +
  geom_ribbon(aes(ymin = f_mean - 2*f_se, ymax = f_mean + 2*f_se, fill = selection), alpha = 0.3) +
  facet_grid(nitrate_grown ~ replicate) +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  labs(x = "Generation", y = "Inbreeding coefficient", tag = "B") +
  theme(legend.position = "none")


# correlation between He and inbreeding
# Calculate median heterozygosity in each population
temp <- pool200 %>%
  filter(sample != "founder population") %>%
  mutate(sample = paste(str_to_lower(nitrate), selection, str_to_lower(rep), sep = "_")) %>%
  group_by(sample) %>%
  summarise(het_mean = mean(hap_het))


p3 <- phen_sum %>%
  filter(generation == 10) %>%
  mutate(sample = paste(str_to_lower(nitrate_grown), selection, str_to_lower(replicate), sep = "_")) %>%
  select(sample, f_mean, selection, nitrate_grown) %>%
  full_join(temp, by = "sample") %>%
  ggplot(aes(f_mean, het_mean, fill = selection)) +
  geom_point(shape = 21, size = 3) +
  labs(x = "Inbreeding coefficient\n(Gen 10)", y = "Heterozygosity", tag = "C") +
  theme(legend.position = "none") +
  scale_y_continuous(breaks = seq(0.8, 0.9, 0.02), limits = c(0.8, 0.9))


# build the figure
p1 / (p2 + p3)
ggsave("./figures/Fig03.pdf", width = 7.5, height = 6)
