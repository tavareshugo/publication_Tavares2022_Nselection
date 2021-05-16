#
# Fig S03
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



#### Fig S03 ####

# SD differential
p1 <- phen_sum %>% 
  ggplot(aes(generation, total_sd_seldif, colour = selection)) +
  geom_line(size = 1, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey") +
  facet_grid(nitrate_grown ~ replicate) +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  scale_y_continuous(breaks = -2:2) +
  labs(x = "Generation", y = "SD differential", tag = "A") +
  scale_color_brewer(palette = "Dark2")

# Median
p2 <- phen_sum %>% 
  filter(selection != "directional") %>% 
  ggplot(aes(generation, total_median, colour = selection)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = total_median - total_mad, ymax = total_median + total_mad, fill = selection), alpha = 0.1,
              colour = NA) + 
  facet_grid(nitrate_grown ~ replicate) +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  labs(x = "Generation", y = "Median branches", tag = "B") +
  scale_color_manual(values = RColorBrewer::brewer.pal(3, name = "Dark2")[c(1, 3)]) + 
  scale_fill_manual(values = RColorBrewer::brewer.pal(3, name = "Dark2")[c(1, 3)])

# SD response
p3 <- var_tests %>% 
  ggplot(aes(generation, sd_dif, colour = replicate)) +
  geom_line(size = 1) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey") +
  facet_grid(nitrate_grown ~ comparison) +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  scale_colour_viridis_d() +
  labs(x = "Generation", y = "SD response\n(selected - control)", tag = "C")


pdf("./figures/FigS03.pdf", width = 7.5, height = 4)
p1 + p2 + plot_layout(ncol = 1)
dev.off()



