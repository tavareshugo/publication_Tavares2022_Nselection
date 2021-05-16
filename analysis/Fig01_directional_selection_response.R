#
# Fig 01
#

#### Setup ####

library(tidyverse)
library(patchwork)
library(rsample)


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
phen_response <- read_rds("./data/processed/phenotypes/phenotypes_response.rds")
phen_sum <- read_rds("./data/processed/phenotypes/phenotypes_summarised.rds")

#### Estimate realised heritability ####

# Make bootstrap resamples stratified per population
resamples <- phen_response %>%
  filter(comparison == "directional") %>%
  select(generation, replicate, nitrate_grown, cumsel, response) %>%
  drop_na() %>%
  mutate(population = interaction(replicate, nitrate_grown)) %>%
  bootstraps(times = 1000, strata = "population")

# Fit model to bootstrap samples
resamples <- resamples %>%
  # fit the model to each sample
  mutate(fit = map(splits, function(i){
    lm(response ~ 0 + replicate:nitrate_grown:cumsel, data = analysis(i)) %>%
      broom::tidy()
  })) %>%
  unnest(fit) %>%
  # tidy the terms
  separate(term, c("replicate", "nitrate", "temp"), sep = ":") %>%
  mutate(replicate = str_remove(replicate, "replicate"),
         nitrate = str_remove(nitrate, "nitrate_grown"))


#### make plots ####

p1 <- phen_sum %>%
  filter(selection != "stabilising") %>%
  # group_by(replicate, nitrate_grown) %>%
  # mutate(total_mean = scale(total_mean, scale = FALSE)) %>%
  # ungroup() %>%
  ggplot(aes(generation, total_median, colour = selection)) +
  geom_line(size = 1) +
  #geom_hline(yintercept = 0, linetype = 2, colour = "grey") +
  geom_ribbon(aes(ymin = total_median - total_mad, ymax = total_median + total_mad, fill = selection), alpha = 0.2,
              colour = NA) +
  facet_grid(nitrate_grown ~ replicate) +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  labs(x = "Generation", y = "Median branches", tag = "A") +
  scale_color_brewer(palette = "Dark2") #+
  #theme(legend.position = c(0, 0.5), legend.justification = c(0, 0.5))

p2 <- phen_response %>%
  filter(comparison == "directional") %>%
  ggplot(aes(cumsel, response, group = replicate, colour = replicate)) +
  geom_point() +
  stat_smooth(method = lm, formula = y ~ 0 + x, se = FALSE) +
  #stat_smooth(method = mgcv::gam, formula = y ~ 0 + x, se = FALSE) +
  facet_grid(nitrate_grown ~ .) +
  geom_line(aes(group = replicate)) +
  labs(x = "Cumulative selection differential",
       y = "Response\n(selected - control)",
       colour = "Replicate",
       tag = "B") +
  theme(legend.position = "none") +
  scale_color_viridis_d()


p3 <- resamples %>%
  ggplot(aes(replicate, estimate)) +
  geom_violin(aes(fill = replicate), scale = "width") +
  facet_grid(nitrate ~ .) +
  scale_y_continuous(limits = c(0, 0.22), breaks = c(0, 0.1, 0.2)) +
  labs(y = "realised heritability", x = "replicate", tag = "C") +
  theme(legend.position = "none") +
  scale_fill_viridis_d() +
  coord_flip()

pdf("./figures/Fig01.pdf", width = 7.5, height = 5)
p1 / (p2 + p3 + plot_layout(widths = c(3/4, 1/4)))
dev.off()
