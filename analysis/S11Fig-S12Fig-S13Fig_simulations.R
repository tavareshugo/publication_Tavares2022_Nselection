#
# Fig S11, S12 and S13
#

library(tidyverse)
library(patchwork)
library(vroom)
library(ggridges)
theme_set(theme_classic() + theme(text = element_text(size = 16)))


# Read data ---------------------------------------------------------------

# heritability
trait_h2 <- vroom("data/processed/simulations/trait_heritabilities.csv")

# inbreeding
trait_sum <- vroom("data/processed/simulations/trait_summarised.csv") %>%
  mutate(selection = factor(selection, levels = c("random", "directional", "stabilising")))

# heterozygosity at selected loci
het_sel <- vroom("data/processed/simulations/het_selected_loci.csv") %>%
  mutate(selection = factor(selection, levels = c("random", "directional", "stabilising")))

# summarised heterozygosity
het_sum <- vroom("data/processed/simulations/het_summarised.csv") %>%
  mutate(selection = factor(selection, levels = c("random", "directional", "stabilising")))

# table of selected loci
selected_loci <- het_sum %>%
  filter(selected) %>%
  distinct(chrom, pos, selected_nloci)

# individual replicates (only reading 10 of them due to memory constraints)
het <- vroom(pipe('cat data/processed/simulations/het_all_loci.csv | grep "seed8[0-9]"'),
             col_names = c("id", 
                           "selected_nloci", "selected_effect",
                           "selected_nalleles", "seed", "selection",
                           "gen", "chrom", "pos", "selected",
                           "het", "het2", "het3", "het4", "het5",
                           "het6", "het7", "het8", "het_quantile"))

# quick plotting filters
filter_nloci <- "35"
filter_nalleles <- 1
filter_effect <- 0.2


# Heritability ------------------------------------------------------------

trait_h2 %>%
  ggplot(aes(factor(selected_nloci), estimate)) +
  annotate(geom = "rect", alpha = 0.3,
           xmin = -Inf, xmax = Inf, ymin = 0.05, ymax = 0.2) +
  geom_violin(aes(fill = factor(selected_effect)), scale = "width") +
  facet_grid(paste0(factor(selected_nalleles), "/19") ~ .) +
  scale_fill_viridis_d(option = "inferno") +
  labs(x = "No. selected loci",
       y = "realised heritability",
       fill = "selected\nallele\neffect")

trait_h2 %>%
  filter(selected_nalleles == filter_nalleles) %>%
  ggplot(aes(factor(selected_nloci), estimate)) +
  annotate(geom = "rect", alpha = 0.3,
           xmin = -Inf, xmax = Inf, ymin = 0.05, ymax = 0.2) +
  geom_violin(aes(fill = factor(selected_effect)), scale = "width") +
  facet_grid(paste0(factor(selected_nalleles), "/19") ~ .) +
  scale_fill_viridis_d(option = "inferno") +
  labs(x = "No. selected loci",
       y = "realised heritability",
       fill = "selected\nallele\neffect")

trait_h2 %>%
  filter(selected_nloci == filter_nloci &
           selected_effect == filter_effect &
           selected_nalleles == filter_nalleles) %>%
  ggplot(aes(estimate)) +
  annotate(geom = "rect", alpha = 0.3,
           xmin = 0.05, xmax = 0.2, ymin = -Inf, ymax = Inf) +
  geom_density(fill = "steelblue") +
  labs(x = "realised heritability")


# Inbreeding --------------------------------------------------------------

trait_sum %>%
  filter(selected_nloci == filter_nloci) %>%
  ggplot(aes(gen, f_mean, fill = selection, group = interaction(selection))) +
  geom_ribbon(stat = "summary", fun.data = "mean_se", fun.args = list(mult = 1.96),
              alpha = 0.5) +
  scale_fill_brewer(palette = "Dark2") +
  facet_grid(selected_effect ~ selected_nalleles)


# Heterozygosity ----------------------------------------------------------

# at selected loci
het_sel %>%
  filter(selected_nloci == filter_nloci &
           selected_nalleles == filter_nalleles &
           selected_effect == filter_effect) %>%
  filter(selection != "stabilising") %>%
  ggplot(aes(selection, het)) +
  geom_violin(aes(fill = selection), scale = "width") +
  scale_fill_brewer(palette = "Dark2") +
  # facet_grid(~ paste(chrom, pos, sep = "-")) +
  labs(x = "Selection",
       y = "He") +
  theme(legend.position = "none")

# summarised
het_sum %>%
  filter(selected_nloci == filter_nloci &
           selected_nalleles == filter_nalleles &
           selected_effect == filter_effect) %>%
  filter(selection != "stabilising") %>%
  ggplot(aes(pos, het_median)) +
  geom_ribbon(aes(ymin = het_q10, ymax = het_q90, fill = selection), alpha = 0.3) +
  geom_line(aes(colour = selection)) +
  geom_point(data = selected_loci %>% filter(selected_nloci == filter_nloci),
             aes(y = 0), fill = "black", size = 3, shape = 24) +
  facet_grid(selected_effect ~ chrom, scales = "free_x", space = "free_x") +
  scale_colour_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Position", y = "Expected heterozygosity", labs = "Selection") +
  theme(axis.text.x = element_blank())


het_sel %>%
  filter(selected_nloci == filter_nloci &
           selected_nalleles == filter_nalleles &
           selected_effect == filter_effect) %>%
  filter(selection != "stabilising") %>%
  left_join(trait_h2) %>%
  ggplot(aes(het, estimate)) +
  geom_point(aes(colour = selection), size = 2) +
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "expected heterozygosity", y = "realised heritability")

# replicates
het %>%
  filter(selection != "stabilising") %>%
  filter(selected_nloci == filter_nloci &
           selected_nalleles == filter_nalleles &
           selected_effect == filter_effect) %>%
  mutate(selection = factor(selection, levels = c("random", "directional", "stabilising"))) %>%
  left_join(trait_h2) %>%
  group_by(seed, selection, chrom) %>%
  # arrange(pos) %>%
  mutate(avg = slider::slide_index_dbl(het, pos, ~ mean(.x), .before = 2, .after = 2)) %>%
  ungroup() %>%
  ggplot(aes(pos, avg)) +
  geom_line(aes(colour = selection), size = 1) +
  geom_point(data = selected_loci %>% filter(selected_nloci == filter_nloci),
             aes(y = 0), colour = "grey", fill = "grey", size = 2, shape = 24) +
  facet_grid(str_remove(seed, "seed20200916") + round(estimate, 2) ~ chrom, scales = "free_x", space = "free_x") +
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "Position", y = "Expected heterozygosity", labs = "Selection",
       title = paste0("loci = ", filter_nloci, "; alleles = ", filter_nalleles, "; effect = ", filter_effect)) +
  theme(axis.text.x = element_blank()) +
  scale_y_continuous(breaks = c(0, 0.5, 1))

het %>%
  filter(selection != "stabilising") %>%
  filter(selected_nloci == filter_nloci &
           selected_nalleles == filter_nalleles &
           selected_effect == filter_effect) %>%
  mutate(selection = factor(selection, levels = c("random", "directional", "stabilising"))) %>%
  left_join(trait_h2) %>%
  mutate(seed = fct_reorder(paste0(seed, selection), het)) %>%
  ggplot(aes(seed, het)) +
  geom_boxplot(aes(fill = selection)) +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "replicate", y = "expected het") +
  theme(axis.text.x = element_blank())



# Illustrative plot -------------------------------------------------------

p1 <- het_sum %>%
  filter(selected_nloci == filter_nloci &
           selected_nalleles == filter_nalleles &
           selected_effect == filter_effect) %>%
  filter(selection != "stabilising") %>%
  ggplot(aes(pos, het_median)) +
  geom_ribbon(aes(ymin = het_q10, ymax = het_q90, fill = selection), alpha = 0.3) +
  geom_line(aes(colour = selection), show.legend = FALSE) +
  geom_point(data = selected_loci %>%
               filter(selected_nloci == filter_nloci),
             aes(y = 0), fill = "black", size = 3, shape = 24) +
  facet_grid( ~ chrom, scales = "free_x", space = "free_x") +
  scale_colour_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  scale_x_continuous(breaks = seq(0, 400, 100)) +
  labs(x = "", y = "Heterozygosity") +
  theme(axis.text.x = element_blank())

p2 <- het_sel %>%
  filter(selected_nloci == filter_nloci &
           selected_nalleles == filter_nalleles &
           selected_effect == filter_effect) %>%
  filter(selection != "stabilising") %>%
  ggplot(aes(selection, het)) +
  geom_violin(aes(fill = selection), scale = "width") +
  scale_fill_brewer(palette = "Dark2") +
  facet_grid(~ paste(chrom, pos, sep = "-")) +
  labs(x = "Selection",
       y = "Heterozygosity",
       subtitle = "at selected loci") +
  theme(legend.position = "none")

p3 <- trait_h2 %>%
  filter(selected_nloci == filter_nloci &
           selected_nalleles == filter_nalleles &
           selected_effect == filter_effect) %>%
  ggplot(aes(estimate)) +
  geom_density(fill = "grey", alpha = 0.5) +
  labs(x = "Heritability")

p4 <- trait_sum %>%
  filter(selected_nloci == filter_nloci &
           selected_nalleles == filter_nalleles &
           selected_effect == filter_effect) %>%
  filter(selection != "stabilising") %>%
  ggplot(aes(gen, f_mean, fill = selection)) +
  geom_ribbon(stat = "summary", fun.data = "mean_se", fun.args = list(mult = 1.96),
              alpha = 0.3) +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Generation", y = "Inbreeding")

p1 /
  ((p2 | p3 | p4) + plot_layout(widths = c(3, 1, 1))) +
  plot_layout(guides = "collect")



# Paper figures ------------------------------------------------------------

# heritability
trait_h2 %>%
  filter(selected_effect %in% c(0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 1)) %>%
  mutate(selected_nloci = factor(as.numeric(selected_nloci))) %>%
  drop_na(selected_nloci) %>%
  ggplot(aes(factor(selected_nloci), estimate)) +
  annotate(geom = "rect", alpha = 0.3,
           xmin = -Inf, xmax = Inf, ymin = 0.05, ymax = 0.2) +
  geom_violin(aes(fill = factor(selected_effect)), scale = "width") +
  facet_grid(paste0(factor(selected_nalleles), "/19") ~ .) +
  scale_fill_viridis_d(option = "inferno") +
  labs(x = "No. trait loci",
       y = "Realised heritability",
       fill = "Allele\nEffect")
  
trait_h2 %>%
  filter(selected_effect %in% c(0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 1),
         selected_nalleles %in% c(1, 4),
         selected_nloci %in% c(1, 5, 10, 30, 60)) %>%
  mutate(selected_nloci = factor(as.numeric(selected_nloci))) %>%
  drop_na(selected_nloci) %>%
  ggplot(aes(estimate, factor(selected_effect))) +
  annotate(geom = "rect", alpha = 0.3,
           ymin = -Inf, ymax = Inf, xmin = 0.05, xmax = 0.2) +
  geom_density_ridges(aes(fill = selected_nloci), 
                      alpha = 0.9) + 
  facet_grid(paste0(selected_nalleles, "/19") ~ .) +
  scale_fill_viridis_d(option = "inferno") +
  labs(x = "Realised Heritability",
       y = "Allelic Effect",
       fill = "No. Adv\nAlleles")
ggsave("./figures/S11Fig.pdf", width = 5, height = 4)

# function to make plot
makePlot <- function(filter_nloci, filter_nalleles, filter_effect){
  p1 <- het_sum %>%
    filter(selected_nloci == filter_nloci &
             selected_nalleles == filter_nalleles &
             selected_effect == filter_effect) %>%
    filter(selection != "stabilising") %>%
    ggplot(aes(pos, het_median)) +
    geom_ribbon(aes(ymin = het_q10, ymax = het_q90, fill = selection), alpha = 0.3) +
    geom_line(aes(colour = selection)) +
    geom_point(data = selected_loci %>% filter(selected_nloci == filter_nloci),
               aes(y = 0), fill = "black", size = 3, shape = 24) +
    facet_grid( ~ chrom, scales = "free_x", space = "free_x") +
    scale_colour_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    labs(x = "Position", y = "He", labs = "Selection",
         subtitle = paste0("Trait loci = ", filter_nloci,
                           "; Start Freq = ", filter_nalleles, "/19",
                           "; Allele effect = ", filter_effect)) +
    theme(axis.text.x = element_blank())

  p2 <- het_sel %>%
    filter(selected_nloci == filter_nloci &
             selected_nalleles == filter_nalleles &
             selected_effect == filter_effect) %>%
    filter(selection != "stabilising") %>%
    filter(chrom == "Chr1") %>%
    ggplot(aes(het)) +
    geom_density(aes(fill = selection), alpha = 0.5) +
    scale_fill_brewer(palette = "Dark2") +
    # facet_grid(~ paste(chrom, pos, sep = "-")) +
    labs(x = "He", y = "") +
    theme(legend.position = "none",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank())

  p3 <- het_sel %>%
    filter(selected_nloci == filter_nloci &
             selected_nalleles == filter_nalleles &
             selected_effect == filter_effect) %>%
    filter(selection != "stabilising") %>%
    group_by(seed, selection, chrom) %>%
    group_by(selection, seed) %>%
    summarise(q1 = sum(het_quantile < 0.01),
              n = n()) %>%
    group_by(selection) %>%
    ungroup() %>%
    ggplot(aes(factor(q1))) +
    geom_bar(aes(fill = selection),
             position = position_dodge2(preserve = "single")) +
    # geom_bar(aes(fill = selection), position = "dodge") +
    scale_fill_brewer(palette = "Dark2") +
    labs(x = "No. loci < Q1", y = "% simulations") +
    theme(legend.position = "none") +
    scale_x_discrete(limits = factor(0:min(10, filter_nloci)),
                     breaks = 0:min(10, filter_nloci))

  (p1 / (p2 | p3))
}

# average across simulations
wrap_plots(
  makePlot(filter_nloci = 10,
           filter_nalleles = 1,
           filter_effect = 0.3) + plot_layout(tag_level = "new"),
  makePlot(filter_nloci = 30,
           filter_nalleles = 1,
           filter_effect = 0.2) + plot_layout(tag_level = "new"),
  ncol = 1
) +
  plot_annotation(tag_levels = c("A", "i"))
ggsave("figures/S12Fig.pdf", width = 7.5, height = 10)

# individual simulation example
filter_nloci <- "30"
filter_nalleles <- 1
filter_effect <- 0.2
het %>%
  filter(selection != "stabilising") %>%
  filter(selected_nloci == filter_nloci &
           selected_nalleles == filter_nalleles &
           selected_effect == filter_effect) %>%
  filter(seed %in% paste0("seed", 80:84)) %>%
  mutate(selection = factor(selection, levels = c("random", "directional", "stabilising"))) %>%
  left_join(trait_h2) %>%
  group_by(seed, selection, chrom) %>%
  # arrange(pos) %>%
  mutate(avg = slider::slide_index_dbl(het, pos, ~ mean(.x), .before = 2, .after = 2)) %>%
  ungroup() %>%
  ggplot(aes(pos, avg)) +
  geom_line(aes(colour = selection), size = 1) +
  geom_point(data = selected_loci %>% filter(selected_nloci == filter_nloci),
             aes(y = 0), fill = "black", size = 2, shape = 24) +
  facet_grid(str_replace(seed, "seed", "sim") + round(estimate, 2) ~ chrom, scales = "free_x", space = "free_x") +
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "Position", y = "He", labs = "Selection") +
  theme(axis.text.x = element_blank()) +
  scale_y_continuous(breaks = c(0, 0.5, 1))
ggsave("figures/S13Fig.pdf", width = 7.5, height = 6)



#### playground ####

filter_nloci <- "30"
filter_nalleles <- 1
filter_effect <- 0.2
het %>%
  filter(selection != "stabilising") %>%
  filter(selected_nloci == filter_nloci &
           selected_nalleles == filter_nalleles &
           selected_effect == filter_effect) %>%
  # filter(seed %in% paste0("seed", 80:85)) %>%
  mutate(selection = factor(selection, levels = c("random", "directional", "stabilising"))) %>% 
  left_join(trait_h2) %>%
  filter(estimate > 0.05) |> 
  group_by(seed, selection, chrom) %>%
  # arrange(pos) %>%
  mutate(avg = slider::slide_index_dbl(het, pos, ~ mean(.x), .before = 2, .after = 2)) %>%
  ungroup() %>%
  ggplot(aes(pos, avg)) +
  geom_line(aes(colour = selection), size = 1) +
  geom_point(data = selected_loci %>% filter(selected_nloci == filter_nloci),
             aes(y = 0), fill = "black", size = 2, shape = 24) +
  facet_grid(str_replace(seed, "seed", "sim") + round(estimate, 2) ~ chrom, scales = "free_x", space = "free_x") +
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "Position", y = "He", labs = "Selection") +
  theme(axis.text.x = element_blank()) +
  scale_y_continuous(breaks = c(0, 0.5, 1))

makePlot(filter_nloci, filter_nalleles, filter_effect) +
   plot_layout(tag_level = "new")


het_sel %>%
  filter(selected_nloci == filter_nloci &
           selected_nalleles == filter_nalleles &
           selected_effect == filter_effect) %>%
  filter(selection != "stabilising") %>%
  filter(selected) %>%
  left_join(trait_h2) %>%
  ggplot(aes(het_quantile, estimate)) +
  geom_point(aes(colour = selection)) +
  geom_hline(yintercept = 0.05) +
  scale_colour_brewer(palette = "Dark2") +
  theme(legend.position = "none")

# heritability broken down
trait_h2 %>%
  mutate(selected_nloci = factor(as.numeric(selected_nloci))) %>%
  drop_na(selected_nloci) %>%
  filter(selected_nalleles == 2) |> 
  ggplot(aes(factor(selected_effect), estimate)) +
  annotate(geom = "rect", alpha = 0.3,
           xmin = -Inf, xmax = Inf, ymin = 0.05, ymax = 0.2) +
  geom_violin(aes(fill = factor(selected_effect)), scale = "width") +
  facet_wrap( ~ selected_nloci) +
  scale_fill_viridis_d(option = "inferno") +
  labs(x = "No. trait loci",
       y = "Realised heritability",
       fill = "Allele\nEffect")

het_sel %>%
  filter(selection != "stabilising") %>%
  filter(selected_nloci == 20 & selected_nalleles == 1 & selected_effect == 0.2) %>%
  group_by(selected_nloci, selected_nalleles, selected_effect,
           selection, seed) %>%
  summarise(q10 = sum(het_quantile < 0.01)) %>%
  ungroup() %>%
  filter(q10 %in% c(0) & selection == "directional") %>%
  left_join(trait_h2) %>%
  filter(estimate > 0.1) %>%
  pull(seed) %>% unique()

het_sel %>%
  filter(selection == "directional") %>%
  mutate(selected_nloci = factor(as.numeric(selected_nloci))) %>%
  drop_na(selected_nloci) %>%
  ggplot(aes(factor(selected_nloci), het_quantile)) +
  geom_violin(aes(fill = factor(selected_effect)), scale = "width") +
  facet_grid(paste0(factor(selected_nalleles), "/19") ~ .) +
  scale_fill_viridis_d(option = "inferno") +
  labs(x = "No. trait loci",
       y = "He quantile",
       fill = "Allele\nEffect")

het %>%
  filter(selection != "stabilising") %>%
  filter(selected_nloci == filter_nloci &
           selected_nalleles == filter_nalleles &
           selected_effect == filter_effect) %>%
  mutate(selection = factor(selection, levels = c("random", "directional", "stabilising"))) %>%
  left_join(trait_h2) %>%
  group_by(seed, selection, chrom) %>%
  # arrange(pos) %>%
  mutate(avg = slider::slide_index_dbl(het, pos, ~ mean(.x), .before = 2, .after = 2)) %>%
  ungroup() %>%
  ggplot(aes(pos, avg)) +
  geom_line(aes(colour = selection), size = 1) +
  geom_point(data = selected_loci %>% filter(selected_nloci == filter_nloci),
             aes(y = 0), colour = "grey", fill = "grey", size = 2, shape = 24) +
  facet_grid(str_remove(seed, "seed20200916") + round(estimate, 2) ~ chrom, scales = "free_x", space = "free_x") +
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "Position", y = "Expected heterozygosity", labs = "Selection",
       title = paste0("loci = ", filter_nloci, "; alleles = ", filter_nalleles, "; effect = ", filter_effect)) +
  theme(axis.text.x = element_blank()) +
  scale_y_continuous(breaks = c(0, 0.5, 1))

het %>%
  filter(selection == "directional") %>%
  filter(selected_nloci == filter_nloci &
           selected_nalleles == filter_nalleles &
           selected_effect == filter_effect) %>%
  mutate(selection = factor(selection, levels = c("random", "directional", "stabilising"))) %>%
  left_join(trait_h2) %>%
  group_by(seed, selection, chrom) %>%
  # arrange(pos) %>%
  mutate(avg = slider::slide_index_dbl(het, pos, ~ mean(.x), .before = 2, .after = 2)) %>%
  mutate(het_quantile2 = ecdf(avg)(avg)) %>%
  group_by(seed, selection) %>%
  mutate(het_quantile3 = ecdf(avg)(avg)) %>%
  ungroup() %>%
  filter(selected) %>%
  select(seed, chrom, pos, het, avg, het_quantile, het_quantile2, het_quantile3) %>%
  arrange(seed, chrom, pos) %>%
  View()

temp <- het %>%
  filter(selected_nloci == filter_nloci &
           selected_nalleles == filter_nalleles &
           selected_effect == filter_effect) %>%
  mutate(selection = factor(selection, levels = c("random", "directional", "stabilising"))) %>%
  left_join(trait_h2) %>%
  group_by(seed, selection, chrom) %>%
  # arrange(pos) %>%
  mutate(avg1 = slider::slide_index_dbl(het, pos, ~ mean(.x), .before = 2, .after = 2),
         avg2 = slider::slide_index_dbl(het, pos, ~ median(.x), .before = 2, .after = 2),
         avg3 = slider::slide_index_dbl(het, pos, ~ min(.x), .before = 2, .after = 2)) %>%
  ungroup()

temp %>%
  select(matches("het"), matches("avg")) %>%
  prcomp(scale = TRUE) %>%
  .$x %>%
  as_tibble() %>%
  bind_cols(temp) %>%
  ggplot(aes(pos, PC1)) +
  geom_line(aes(colour = selection), size = 1) +
  geom_point(data = selected_loci %>% filter(selected_nloci == filter_nloci),
             aes(y = 0), colour = "grey", fill = "grey", size = 2, shape = 24) +
  facet_grid(str_remove(seed, "seed20200916") + round(estimate, 2) ~ chrom, scales = "free_x", space = "free_x") +
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "Position", y = "Expected heterozygosity", labs = "Selection",
       title = paste0("loci = ", filter_nloci, "; alleles = ", filter_nalleles, "; effect = ", filter_effect)) +
  theme(axis.text.x = element_blank())

