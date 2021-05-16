library(slider)
library(tidyverse)
library(patchwork)
theme_set(theme_classic() + theme(text = element_text(size = 16)))

library(reticulate)
use_condaenv("simupop")
source_python("workflow/scripts/simulations/single_replicate_simulation.py")
rm(list = ls())


sim_het <- read_csv("heterozygosity.csv")

sim_het <- sim_het %>%
  filter(selection != "stabilising") %>%
  separate(locus, into = c("chrom", "pos"), sep = "-",
           remove = FALSE, convert = TRUE) %>%
  mutate(selection = factor(selection, levels = c("random", "directional", "stabilising")))

selected_loci <- sim_het %>% filter(selected) %>% distinct(selection, chrom, pos)

p1 <- sim_het %>%
  filter(gen %in% c(10)) %>%
  rename(avg = het) %>%
  # group_by(selection, gen, chrom) %>%
  # mutate(avg = slide_index_dbl(het, pos, ~ mean(.x), .before = 5, .after = 5)) %>%
  # ungroup() %>%
  ggplot(aes(pos, avg)) +
  geom_point(data = selected_loci, aes(y = 1),
             fill = "black", pch = 25, size = 4) +
  #geom_line(aes(y = het), colour = "lightgrey") +
  geom_line(aes(colour = selection)) +
  facet_grid( ~ chrom, scale = "free_x", space = "free_x") +
  scale_colour_brewer(palette = "Dark2") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 2000, 500))

p2 <- sim_het %>%
  filter(gen %in% c(0, 10)) %>%
  ggplot(aes(factor(gen), het)) +
  geom_boxplot(aes(fill = selection)) +
  geom_hline(yintercept = 1 - sum(19 * 1/19^2)) +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Generation", y = "He")




# Inbreeding --------------------------------------------------------------
# not loading the whole package due to conflicting functions
inbreeding <- function(ind, mother, father){
  d <- data.frame(ind, mother, father)
  pedigree::calcInbreeding(d)
}

# read pedigree
sim_ped <- read_csv("pedigree.csv") %>%
  mutate(selection = factor(selection, levels = c("random", "directional", "stabilising")),
         selected = fitness == 1)

# sim_ped %>%
#   ggplot(aes(trait, rank)) +
#   geom_point(aes(colour = factor(gen), shape = factor(fitness))) +
#   scale_colour_viridis_d() +
#   scale_shape_manual(values = c(20, 1)) +
#   facet_grid(~ selection)
#
# # check trait value
# sim_ped %>%
#   filter(gen %in% c(0, 5, 10)) %>%
#   ggplot(aes(factor(nalleles), trait)) +
#   geom_hline(yintercept = 0, lty = 2) +
#   geom_boxplot(aes(fill = selection)) +
#   facet_grid(gen ~ .)
#
# sim_ped %>%
#   ggplot(aes(factor(nalleles), trait)) +
#   geom_hline(yintercept = 0, lty = 2) +
#   geom_boxplot()


# inbreeding
p3 <- sim_ped %>%
  group_by(selection) %>%
  mutate(f = inbreeding(ind_id, mother_id, father_id)) %>%
  ggplot(aes(factor(gen), f)) +
  geom_ribbon(stat = "summary", fun.data = "mean_se",
              aes(group = selection, fill = selection),
              alpha = 0.5) +
  # geom_line(stat = "summary", fun.data = "mean_se", aes(group = selection)) +
  scale_fill_brewer(palette = "Dark2")

# p1 /
#   ((p2 | p3) + plot_layout(widths = c(2, 1))) +
#   plot_annotation(title = "Simulated data") +
#   plot_layout(guides = "collect")



# Trait response ----------------------------------------------------------

# 1 row per generation, nitrate, selection, replicate, population
# Calculate summary metrics and also selection differentials
sim_ped_summary <- sim_ped %>%
  group_by(gen, selection) %>%
  summarise_at(vars(trait),
               list(~ mean(., na.rm = TRUE),
                    ~ sd(., na.rm = TRUE),
                    se = ~ sd(., na.rm = TRUE)/sqrt(sum(!is.na(.))),
                    ~ median(., na.rm = TRUE),
                    ~ mad(., constant = 1, na.rm = TRUE),
                    q25 = ~ quantile(., na.rm =  TRUE, 0.25),
                    q75 = ~ quantile(., na.rm =  TRUE, 0.75),
                    iqr = ~ IQR(., na.rm = TRUE),
                    n = ~ sum(!is.na(.)),
                    mean_sel = ~ mean(.[selected], na.rm = TRUE),
                    sd_sel = ~ sd(.[selected], na.rm = TRUE),
                    se_sel = ~ sd(.[selected], na.rm = TRUE)/sqrt(sum(!is.na(.[selected]))),
                    median_sel = ~ median(.[selected], na.rm = TRUE),
                    mad_sel = ~ mad(.[selected], constant = 1, na.rm = TRUE),
                    iqr_sel = ~ IQR(.[selected], na.rm = TRUE),
                    n_sel = ~ sum(!is.na(.[selected])))) %>%
  group_by(selection) %>%
  arrange(gen) %>%
  mutate(seldif = mean_sel - mean,
         cumsel = lag(cumsum(seldif)),
         sd_seldif = sd_sel - sd) %>%
  ungroup()

phen_response <- sim_ped_summary %>%
  group_by(gen) %>%
  summarise(cumsel_directional = cumsel[selection == "directional"],
            response_directional = mean[selection == "directional"] - mean[selection == "random"],
            cumsel_stabilising = cumsel[selection == "stabilising"],
            response_stabilising = mean[selection == "stabilising"] - mean[selection == "random"]) %>%
  ungroup() %>%
  gather(key, value, cumsel_directional:response_stabilising) %>%
  separate(key, c("metric", "comparison")) %>%
  spread(metric, value)

# Make bootstrap resamples stratified per population
library(rsample)
resamples <- phen_response %>%
  filter(comparison == "directional") %>%
  select(gen, cumsel, response) %>%
  drop_na() %>%
  bootstraps(times = 1000)

# Fit model to bootstrap samples
resamples <- resamples %>%
  # fit the model to each sample
  mutate(fit = map(splits, function(i){
    lm(response ~ 0 + cumsel, data = analysis(i)) %>%
      broom::tidy()
  })) %>%
  unnest(fit)

p4 <- phen_response %>%
  filter(comparison == "directional") %>%
  ggplot(aes(cumsel, response)) +
  geom_point() +
  stat_smooth(method = lm, formula = y ~ 0 + x, se = FALSE) +
  geom_line() +
  labs(x = "Cumulative selection differential",
       y = "Response\n(selected - control)")

p5 <- resamples %>%
  ggplot(aes(estimate)) +
  geom_density() +
  labs(x = "realised heritability") +
  coord_cartesian(xlim = c(0, 0.5))

p1 /
  ((p2 | p3) + plot_layout(widths = c(2, 1))) /
  ((p4 | p5) + plot_layout(widths = c(2, 1))) +
  plot_annotation(title = "Simulated data") +
  plot_layout(guides = "collect")


# expectedHomozygosity <- function(freq, pool_alleles = 1){
#
#   # check input
#   if(!is.numeric(freq)) stop("freq must be numeric vector.")
#   if(pool_alleles < 1 | pool_alleles >= length(freq)){
#     warning("pool_alleles is smaller or equal to number of alleles. Assuming hom = 1")
#     return(1)
#   } else{
#     # sort frequencies
#     freq <- sort(freq, decreasing = TRUE)
#
#     # get alleles to pool together
#     freq1 <- sum(freq[1:pool_alleles])
#     freq2 <- freq[-(1:pool_alleles)]
#
#     # expected homozygosity
#     hom <- freq1^2 + sum(freq2^2)
#
#     return(hom)
#   }
# }
#
# sim_freq <- read_csv("frequency.csv")
# sim_freq <- sim_freq %>%
#   filter(selection != "stabilising") %>%
#   separate(locus, into = c("chrom", "pos"), sep = "-",
#            remove = FALSE, convert = TRUE) %>%
#   mutate(selection = factor(selection, levels = c("random", "directional", "stabilising"))) %>%
#   group_by(selection, generation, locus) %>%
#   mutate(het = 1 - expectedHomozygosity(freq, 1),
#          het2 = 1 - expectedHomozygosity(freq, 2),
#          het3 = 1 - expectedHomozygosity(freq, 3),
#          het4 = 1 - expectedHomozygosity(freq, 4)) %>%
#   ungroup()
#
# sim_het %>%
#   filter(gen %in% c(10)) %>%
#   rename(avg = het) %>%
#   # group_by(selection, gen, chrom) %>%
#   # mutate(avg = slide_index_dbl(het, pos, ~ mean(.x), .before = 5, .after = 5)) %>%
#   # ungroup() %>%
#   ggplot(aes(pos, avg)) +
#   geom_point(data = selected_loci, aes(y = 1),
#              fill = "black", pch = 25, size = 4) +
#   #geom_line(aes(y = het), colour = "lightgrey") +
#   geom_line(aes(colour = selection)) +
#   facet_grid( ~ chrom, scale = "free_x", space = "free_x") +
#   scale_colour_brewer(palette = "Dark2") +
#   scale_y_continuous(limits = c(0, 1)) +
#   scale_x_continuous(breaks = seq(0, 2000, 500))
