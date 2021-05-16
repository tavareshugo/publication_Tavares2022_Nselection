#
# Fig 04
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
source("./analysis/data_processing/read_phenotypes.R")

# heterozygosity scans at 200Kb windows
pool200 <- read_csv("./data_processed/poolseq/poolseq_window200kb.csv",
                    col_types = "ccccciiiiddddddddccddd") %>%
  mutate(selection = factor(selection, levels = c("random", "directional", "stabilising")))

# read starting population scans at 200kb windows
startpop200 <- read_csv("./data_processed/startpop/startpop_window200kb.csv",
                        col_types = "iiiiidididdddddddddddddddddd")

# put these two together
pool200 <- startpop200 %>%
  mutate(sample = "founder population") %>%
  select(sample, chrom, start = win_start, end = win_end, hap_het = hap_het_mean) %>%
  bind_rows(pool200)

#### make graphs ####

# boxplot of heterozygosity
p1 <- pool200 %>%
  group_by(sample) %>%
  mutate(med = median(hap_het)) %>%
  ungroup() %>%
  mutate(sample = reorder(sample, med)) %>%
  ggplot(aes(sample, hap_het, fill = selection, group = sample)) +
  geom_boxplot(position = position_dodge2(preserve = "total")) +
  geom_point(stat = "summary", fun.y = "mean") +
  geom_hline(yintercept = 1 - (19 * ((1/19)^2)),
             linetype = "dashed", colour = "grey48") +
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
pdf("./figures/Fig04.pdf", width = 7.5, height = 6)
p1 / (p2 + p3)
dev.off()
