#
# Fig 02
#

#### Setup ####

library(tidyverse)
library(patchwork)
library(here)
library(lme4)
library(emmeans)

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

# Read data
gen10_plast <- read_rds("./data/processed/phenotypes/generation10_plasticity.rds")
phen_sum <- read_rds("./data/processed/phenotypes/phenotypes_summarised.rds")
phen_filter <- read_rds("./data/processed/phenotypes/phenotypes_individual.rds")

# source custom function
source("./analysis/functions/wilcoxonRankSum.R")


#### Calculate response flowering #####

# Calculate difference between stabilising of selected and random populations
bolt_response <- phen_sum %>%
  group_by(generation, replicate, nitrate_grown) %>%
  summarise(cumsel_directional = bolt_cumsel[selection == "directional"],
            response_directional = bolt_mean[selection == "directional"] - bolt_mean[selection == "random"],
            cumsel_stabilising = bolt_cumsel[selection == "stabilising"],
            response_stabilising = bolt_mean[selection == "stabilising"] - bolt_mean[selection == "random"]) %>%
  ungroup() %>%
  gather(key, value, cumsel_directional:response_stabilising) %>%
  separate(key, c("metric", "comparison")) %>%
  spread(metric, value)


# Run non-parametric comparison and obtain effect sizes
## Comparing "directional" with "random"
directional_test <- phen_filter %>%
  mutate(id = 1:n()) %>%
  rename(bolt = sow_to_bolt) %>%
  select(id, generation, replicate, nitrate_grown, bolt, selection) %>%
  filter(complete.cases(.)) %>%
  spread(selection, bolt) %>%
  group_by(generation, nitrate_grown, replicate) %>%
  do(wilcoxonRankSum(.$directional, .$random, conf.int = TRUE)) %>%
  mutate(comparison = "directional") %>%
  ungroup()

## comparing "stabilising" with "random"
stabilising_test <- phen_filter %>%
  mutate(id = 1:n()) %>%
  rename(bolt = sow_to_bolt) %>%
  select(id, generation, replicate, nitrate_grown, bolt, selection) %>%
  filter(complete.cases(.)) %>%
  spread(selection, bolt) %>%
  group_by(generation, nitrate_grown, replicate) %>%
  do(wilcoxonRankSum(.$stabilising, .$random, conf.int = TRUE)) %>%
  mutate(comparison = "stabilising") %>%
  ungroup()

## bind the results in a single table
all_tests <- bind_rows(directional_test, stabilising_test) %>%
  ungroup() %>%
  mutate(p.adjust = p.adjust(p.value, method = "fdr"))

## Join with main response table
bolt_response <- full_join(bolt_response, all_tests,
                           by = c("generation", "replicate", "nitrate_grown", "comparison"))

rm(directional_test, stabilising_test, all_tests)


#### Calculate plasticity ####

phen <- read_csv("./data/raw/phenotypes/compiled_phenotypes_data_clean.csv",
                 col_types = "iccccDiiinninccccccclnDiiiDnnncccllc")

# Change selection names to be more informative
phen <- phen %>%
  mutate(selection = case_when(selection == "random" ~ "random",
                               selection == "most" ~ "directional",
                               selection == "average" ~ "stabilising",
                               TRUE ~ as.character(NA))) %>%
  mutate(selection = factor(selection,
                            levels = c("random", "directional", "stabilising")),
         nitrate_grown = factor(nitrate_grown, levels = c("LN", "HN")))


gen10 <- phen %>%
  filter(generation == 10 & !duplicated_line & !population_discordant) %>%
  select(population, selection, nitrate_selected,
         nitrate_grown, replicate, family, total)


# fit model separately for each set of populations
# (I couldn't figure out how to implement a single model with right levels of nesting)
contrast_output <- list()
for(i in c("LN", "HN")){
  for(j in c("A", "B", "C")){
    contrast_output[[paste(i, j, sep = "_")]] <- gen10 %>%
      # subset relevant part of the data
      filter(nitrate_selected == i & replicate == j) %>%
      # fit model - accounts for nested design
      lmer(total ~ selection*nitrate_grown + (nitrate_grown|family),
           data = .) %>%
      emmeans( ~ nitrate_grown*selection) %>%
      contrast(list(`Directional vs Random` = c(1, -1, -1, 1, 0, 0),
                    `Stabilising vs Random` = c(1, -1, 0, 0, -1, 1)),
               adjust = "none") %>%
      as_tibble()
  }
}

contrast_output %>%
  bind_rows(.id = "set") %>%
  separate(set, c("nitrate_grown", "replicate"), sep = "_") %>%
  mutate(padj = p.adjust(p.value, method = "fdr")) %>%
  arrange(nitrate_grown, contrast, replicate) %>%
  mutate(print = paste0("replicate ", replicate, ": ", round(estimate, 2), " +/- ", round(SE, 2), ", p = ", round(padj, 4), "; ")) %>%
  select(nitrate_grown, contrast, print)


#### Paper figure ####

# Plasticity response
p1 <- gen10_plast %>%
  ggplot(aes(selection, total_mean_D, colour = selection)) +
  geom_hline(yintercept = 0, colour = "grey48", linetype = "dashed") +
  # ggbeeswarm::geom_quasirandom(colour = "grey", size = 3) +
  geom_dotplot(binaxis = "y", method = "histodot",
               stackdir = "center", fill = "grey", colour = "grey", binwidth = 1, dotsize = 0.8) +
  geom_pointrange(stat = "summary", fun.data = "mean_cl_boot") +
  facet_grid(nitrate_selected ~ replicate) +
  labs(x = "Selection regime", y = "Plasticity\n(HN - LN)", tag = "A") +
  scale_y_continuous(breaks = seq(-4, 6, 2)) +
  theme(legend.position = "none")

# Flowering response
p2 <- bolt_response %>%
  filter(comparison == "directional") %>%
  ggplot(aes(generation, vda, colour = replicate)) +
  geom_line(size = 1) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey") +
  facet_grid(nitrate_grown ~ .) +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  labs(x = "Generation", y = "Flowering response\n(A statistic)", tag = "B") +
  scale_color_viridis_d()

p1 + p2 + plot_layout(ncol = 1)
ggsave("./figures/Fig02.pdf", width = 7.5, height = 5)
