#
# Fig S06
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

source("./analysis/functions/corBoot.R")

# Read data
magics <- read_csv("./data/raw/phenotypes/compiled_phenotypes_data_clean.csv",
                 col_types = "iccccDiiinninccccccclnDiiiDnnncccllc")
phen <- readRDS("data/processed/phenotypes/phenotypes_individual.rds")
phen_sum <- read_rds("./data/processed/phenotypes/phenotypes_summarised.rds")


#### summarise MAGIC data ####

# Summarise data in MAGIC lines
magic_sum <- magics %>%
  # get MAGIC line founders (generation "-1")
  filter(generation == -1 & line %in% c(mother, father)) %>%
  rename(bolt = sow_to_bolt) %>%
  group_by(line, nitrate_grown) %>%
  summarise(across(c(bolt, height, total),
                   list(
                     mean = ~mean(., na.rm = TRUE),
                     sd = ~sd(., na.rm = TRUE),
                     median = ~median(., na.rm = TRUE),
                     mad = ~ mad(., constant = 1),
                     n = ~ sum(!is.na(.))
                   )
  )
  ) %>%
  ungroup()

# Trait correlations in the magic lines
magic_cor <- magic_sum %>%
  group_by(nitrate_grown) %>%
  nest() %>%
  mutate(bolt_total = map(data, ~ corBoot(., formula = "~ bolt_median + total_median", ncpus = 2,
                                          method = "spearman", use = "complete.obs"))) %>%
  unnest(bolt_total, .sep = "_")

# Trait correlations across the generations
# bootstrap CI take a while to run
sel_cor <- phen %>%
  rename(bolt = sow_to_bolt) %>%
  group_by(generation, replicate, nitrate_grown, selection) %>%
  nest() %>%
  mutate(bolt_total = map(data, ~ corBoot(., formula = "~ bolt + total", ncpus = 2,
                                          method = "spearman", use = "complete.obs")))


#### Paper figure ####

# MAGIC
p1 <- magic_sum %>%
  ggplot(aes(total_median, bolt_median)) +
  geom_point(alpha = 0.5) +
  #geom_smooth(se = FALSE) +
  geom_text(data = magic_cor, x = 9, y = 80, hjust = 1, size = 3,
            aes(label = paste0("rho = ", round(bolt_total_statistic, 2),
                               " [", round(bolt_total_conflo, 2), ", ",
                               round(bolt_total_confhi, 2), "]"))) +
  facet_grid(. ~ nitrate_grown) +
  labs(x = "Total branches", y = "Days to flowering", tag = "A")

# Correlation across generations
p2 <- sel_cor %>%
  filter(selection != "stabilising") %>%
  unnest(bolt_total) %>%
  ggplot(aes(generation, statistic)) +
  geom_line(aes(colour = selection)) +
  geom_ribbon(aes(ymin = conflo, ymax = confhi, fill = selection), alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey48") +
  scale_x_continuous(breaks = 0:10) +
  facet_grid(nitrate_grown ~ replicate) +
  labs(x = "Generation", y = "Spearman correlation\n(flowering vs. branches)",
       tag = "B")


# Selection intensity
p3 <- phen_sum %>%
  filter(selection != "stabilising") %>%
  ggplot(aes(generation, bolt_seldif/bolt_sd)) +
  geom_line(aes(colour = selection), size = 1) +
  facet_grid(nitrate_grown ~ replicate) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey48") +
  scale_x_continuous(breaks = 0:10) +
  labs(x = "Generation", y = "Selection intensity", tag = "C")

# Response to selection
p4 <- phen_sum %>%
  filter(selection != "stabilising") %>%
  ggplot(aes(generation, bolt_median, colour = selection)) +
  geom_line() +
  geom_ribbon(aes(ymin = bolt_median - bolt_mad, ymax = bolt_median + bolt_mad,
                  colour = NULL, fill = selection),
              alpha = 0.3) +
  facet_grid(nitrate_grown ~ replicate) +
  labs(x = "Generation", y = "Days to flowering", tag = "D") +
  scale_x_continuous(breaks = 0:10) +
  theme(legend.position = "none")

p1 + p2 + p3 + p4 + plot_layout(ncol = 1)
ggsave("./figures/S06Fig.pdf", width = 7.5, height = 8.70)
