#
# Fig S01
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

# Read data
phen_sum <- readRDS("./data/processed/phenotypes/phenotypes_summarised.rds")


#### Figure ####

# Selection differential
p1 <- phen_sum %>% 
  filter(selection != "stabilising") %>% 
  ggplot(aes(generation, total_seldif)) +
  geom_line(aes(colour = selection), size = 1) +
  facet_grid(nitrate_grown ~ replicate) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey48") +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  labs(x = "Generation", y = "Selection differential", tag = "B")

# Relationship with variance
p2 <- phen_sum %>% 
  filter(selection == "directional") %>% 
  ggplot(aes(total_sd, total_seldif, colour = replicate)) +
  geom_point(size = 1) +
  facet_grid(nitrate_grown ~ .) +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  labs(x = "Standard deviation", y = "Selection differential", tag = "C") +
  scale_colour_viridis_d() +
  theme(legend.position = c(1, 0), legend.justification = c(1, 0))

# Selection intensity
p3 <- phen_sum %>% 
  filter(selection != "stabilising") %>% 
  ggplot(aes(generation, total_seldif/total_sd)) +
  geom_line(aes(colour = selection), size = 1) +
  facet_grid(nitrate_grown ~ replicate) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey48") +
  scale_x_continuous(breaks = 0:10) +
  labs(x = "Generation", y = "Selection intensity", title = "D") +
  theme(legend.position = "none")

p1 / (p2 + p3 + plot_layout(widths = c(1/4, 3/4)))
ggsave("./figures/S01Fig.pdf", width = 7.5, height = 5)
