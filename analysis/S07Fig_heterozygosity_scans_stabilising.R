#
# Fig S07
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

# This reads the pool-scans, centromere location and candidate sweep regions
source("./analysis/data_processing/read_poolseq.R")

# filter candidate sweeps to retain only ones of interest for this figure
candidate_sweeps <- candidate_sweeps %>% 
  filter(selection != "directional" & statistic == "H1")

#### make graphs ####

p1 <- pool200 %>% 
  filter(selection != "directional") %>% 
  ggplot() +
  geom_line(aes(pos/1e6, hap_het, colour = selection)) +
  geom_hline(yintercept = threshold$hap_hom1_q1, linetype = 2, colour = "grey48") +
  geom_rect(data = centromeres,
            aes(xmin = start/1e6 - 0.5, xmax = end/1e6 + 0.5, ymin = -Inf, ymax = Inf), fill = "grey", alpha = 0.5) +
  geom_segment(data = candidate_sweeps,
               aes(x = start/1e6, xend = end/1e6, y = 0.3, yend = 0.3), size = 3) +
  facet_grid(nitrate + rep ~ chrom, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = seq(0, 30, 10)) +
  labs(x = "Mb", y = "Heterozygosity", colour = "selection: ") +
  theme(legend.position = "top") +
  scale_color_manual(values = RColorBrewer::brewer.pal(3, name = "Dark2")[c(1, 3)])

# save graph
p1
ggsave("./figures/S07Fig.pdf", width = 7.5, height = 8.5)
