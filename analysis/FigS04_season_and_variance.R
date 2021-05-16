#
# Fig S04
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



#### Fig S04 ####

# An example of high CV generation
temp <- phen_filter %>% 
  filter(generation == 9 & !is.na(total)) %>% 
  group_by(selection, replicate, nitrate_grown, sow_month) %>% 
  summarise(cv = sd(total)/mean(total))

p1 <- phen_filter %>% 
  filter(generation == 9 & !is.na(total)) %>% 
  ggplot(aes(total)) +
  geom_histogram(binwidth = 1,
                 aes(fill = selection), position = "dodge") +
  facet_grid(nitrate_grown ~ paste0(replicate, "\n(", month.name[sow_month], ")")) +
  geom_text(data = temp, x = 12, y = rep(c(80, 100, 120), 6), 
            aes(label = round(cv, 2), colour = selection), 
            hjust = 0, size = 3) +
  labs(x = "Branch number", tag = "A") +
  theme(legend.position = "none")
theme(legend.position = "none")


# Distribution of shoot branching grouped by sowing month
p2 <- phen_filter %>% 
  filter(!is.na(sow_month) & !is.na(total)) %>% 
  group_by(nitrate_grown, selection, sow_month) %>% 
  mutate(cv = sd(total, na.rm = TRUE)/mean(total, na.rm = TRUE)) %>% 
  ungroup() %>% 
  ggplot(aes(x = total, y = factor(sow_month), fill = cv)) +
  ggridges::geom_density_ridges_gradient(stat = "binline", binwidth = 1) + 
  facet_grid(nitrate_grown ~ selection) +
  scale_y_discrete(labels = month.abb) +
  scale_fill_continuous(trans = "log2") +
  labs(y = "Sowing month", x = "Number of branches", tag = "B")


pdf("./figures/FigS04.pdf", width = 7.5, height = 6)
p1 + p2 + plot_layout(ncol = 1, heights = c(1, 2))
dev.off()

