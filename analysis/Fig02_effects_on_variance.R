#
# Fig 02
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
var_tests <- read_rds("./data/processed/phenotypes/variance_tests.rds")
phen_sum <- read_rds("./data/processed/phenotypes/phenotypes_summarised.rds")
phen_filter <- read_rds("./data/processed/phenotypes/phenotypes_individual.rds")

##### Fig 02 ####

# CV response
p1 <- var_tests %>% 
  ggplot(aes(generation, cv_dif, colour = replicate)) +
  geom_line(size = 1) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey") +
  facet_grid(nitrate_grown ~ comparison) +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  scale_colour_viridis_d() +
  labs(x = "Generation", y = "CV response\n(selected - control)", tag = "A")

# CV across generations
p2 <- phen_sum %>% 
  ggplot(aes(generation, total_sd/total_mean, colour = selection)) +
  geom_line(size = 1) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey") +
  facet_grid(nitrate_grown ~ replicate) +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  scale_y_continuous(trans = "log2", breaks = c(0.3, 0.5, 1, 2, 3)) +
  labs(x = "Generation", y = "CV", tag = "B")

# variance and season
p3 <- phen_filter %>% 
  filter(!is.na(sow_month) & !is.na(total)) %>% 
  group_by(generation, selection, nitrate_grown, replicate, sow_month, population) %>% 
  summarise(total_cv = sd(total)/mean(total), n = n()) %>% 
  filter(n >= 10) %>% 
  ggplot(aes(as.numeric(sow_month), total_cv, colour = selection)) +
  geom_point(position = position_dodge(0.4)) +
  geom_smooth(se = FALSE) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey") +
  scale_x_continuous(breaks = 1:12, labels = month.abb[1:12]) + 
  scale_y_continuous(trans = "log2", breaks = c(0.3, 0.5, 1, 2, 3)) +
  facet_grid(nitrate_grown ~ .) +
  labs(x = "Sowing month", y = "CV", tag = "C") +
  


pdf("./figures/Fig02.pdf", width = 7.5, height = 7)
p1 + p2 + p3 + plot_layout(ncol = 1)
dev.off()

