#
# Fig 01
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
source("./scripts/R/data_processing/read_phenotypes.R")


#### make plots ####

p1 <- phen_response %>% 
  filter(comparison == "directional") %>% 
  ggplot(aes(generation, vda, colour = replicate)) +
  geom_line(size = 1) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey") + 
  facet_grid(nitrate_grown ~ .) +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  labs(x = "Generation", y = "A") +
  scale_color_viridis_d() + 
  theme(legend.position = c(0, 1), legend.justification = c(0, 1))

pdf("./figures/FigS02.pdf", width = 7.5, height = 4)
p1
dev.off()
