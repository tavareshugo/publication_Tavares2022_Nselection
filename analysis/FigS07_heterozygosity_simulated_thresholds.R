library(tidyverse)
library(ggridges)
library(patchwork)

theme_set(theme_classic())


# Read pedigree -----------------------------------------------------------

# get experimental pedigree
ped <- read_csv("./data/raw/phenotypes/compiled_phenotypes_data_clean.csv",
                col_types = "iccccDiiinninccccccclnDiiiDnnncccllc")

# remove few individuals without complete pedigree information
ped <- ped %>%
  distinct(generation, line, mother, father)

# there are some founders with missing phenotype data, so we add them separately
miss <- ped %>%
  filter(!(father %in% line & mother %in% line) & generation != -1) %>%
  select(mother, father) %>%
  pivot_longer(c(mother, father), values_to = "line") %>%
  distinct(line) %>%
  mutate(generation = -1L, mother = NA, father = NA) %>%
  select(generation, line, mother, father)

ped <- bind_rows(miss, ped) %>%
  distinct(generation, line, mother, father) %>%
  arrange(generation, line, mother, father)

rm(miss)



# Read other data ---------------------------------------------------------

# empirical heterozygosity (to compare with simulations)
het <- read_csv("data/processed/poolseq/haplotype_diversity_200kb.csv") %>%
  mutate(sample = str_remove(sample, "_reseq"),
         het = 1 - hap_hom1)

# centromeres (to annotate plot)
centromeres <- read_csv("data/external/tair9_centromeres.csv")


# Simulation --------------------------------------------------------------

# initiate a matrix of genotypes
n_sim_loci <- 1000 # number of loci to simulate
genotypes <- matrix(nrow = nrow(ped), ncol = n_sim_loci*2)
rownames(genotypes) <- ped$line
colnames(genotypes) <- paste0("locus", rep(1:n_sim_loci, each = 2), ".", c(1, 2))

# simulate some parental genotypes
sim_genos <- function(nloci = n_sim_loci){
  alleles <- c("bur", "can", "col", "ct", "edi", "hi", "kn", "ler", "mt", "no", "oy", "po", "rsch", "sf", "tsu", "wil", "ws", "wu", "zu")

  out <- sample(alleles, nloci, replace = TRUE)
  out <- rep(out, each = 2) # homozygous
  return(out)
}

set.seed(1596467939) # for reproducibility
system.time(
  for(i in 1:nrow(ped)){
    # get IDs of the current individual and its parents
    id <- ped$line[i]
    mother <- ped$mother[i]
    father <- ped$father[i]

    cat(i, id, mother, father, "\n") # print message for tracking progress

    # sometimes parents are repeated so we should only sample them once
    if(any(is.na(genotypes[id, ]))){
      if(any(!is.na(genotypes[id, ]))) stop("Some NA's!?")

      # simulate genotype if it's the founder generation
      if(ped$generation[i] == -1){
        genotypes[id, ] <- sim_genos()

      } else { # otherwise sample from the parents
        maternal_alleles <- paste0("locus", 1:n_sim_loci, ".",
                                   sample(c(1, 2), n_sim_loci, replace = TRUE))
        paternal_alleles <- paste0("locus", 1:n_sim_loci, ".",
                                   sample(c(1, 2), n_sim_loci, replace = TRUE))
        # nice trick: https://stackoverflow.com/a/25961969/5023162
        genotypes[id, ] <- c(rbind(genotypes[mother, maternal_alleles],
                                   genotypes[father, paternal_alleles]))
      }
    }
  }
)

# clear workspace
rm(i, id, mother, father, maternal_alleles, paternal_alleles)



# Calculate simulated heterozygosity --------------------------------------

# read processed phenotype data to fetch generation 10 individuals in each population
gen10_idx <- read_rds("data/processed/phenotypes/phenotypes_individual.rds") %>%
  select(generation, line, nitrate_grown, selection, replicate,
         line, mother, father) %>%
  filter(generation == 10) %>%
  mutate(population = paste(nitrate_grown, selection, replicate, sep = "-")) %>%
  with(split(line, population))

# calculate allele frequencies in generation 10
sim_freqs <- lapply(gen10_idx, function(ids){
  alleles <- c("bur", "can", "col", "ct", "edi", "hi", "kn", "ler", "mt", "no", "oy", "po", "rsch", "sf", "tsu", "wil", "ws", "wu", "zu")

  # sample from the genotypes to account for PoolSeq sampling step
  # assuming depth of coverage 100x
  sampled_ids <- sample(ids, 100, replace = TRUE)

  out <- matrix(nrow = 19, ncol = ncol(genotypes))
  rownames(out) <- alleles
  colnames(out) <- colnames(genotypes)

  for(i in alleles){
    out[i, ] <- colSums(genotypes[sampled_ids, ] == i)/length(sampled_ids)
  }
  return(out)
})

# tidy things up
sim_freqs <- map_dfr(sim_freqs, function(d){
  d %>%
    as_tibble(rownames = "allele") %>%
    pivot_longer(-"allele", names_to = "locus", values_to = "count")
}, .id = "population")

# summarise the two copies of each locus as one
sim_freqs <- sim_freqs %>%
  mutate(locus = str_remove(locus, "locus")) %>%
  separate(locus, c("locus", "copy"), convert = TRUE) %>%
  group_by(population, allele, locus) %>%
  summarise(count = sum(count)/2) %>%
  ungroup()

# calculate simulated heterozygosity
sim_het <- sim_freqs %>%
  group_by(locus, population) %>%
  summarise(het = 1 - sum(count^2),
            tot = sum(count)) %>%
  ungroup() %>%
  separate(population, c("nitrate", "selection", "rep"), remove = FALSE)


# Visualise -----------------------------------------------------

# distribution of simulated heterozygosity
p1 <- sim_het %>%
  mutate(population = str_replace_all(population, "-", "\n")) %>%
  mutate(population = fct_reorder(population, het)) %>%
  mutate(selection = factor(selection, levels = c("random", "directional", "stabilising"))) %>%
  ggplot(aes(population, het)) +
  geom_boxplot(aes(group = population, fill = selection)) +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Population", y = "Heterozygosity\n(simulated)")

# compare with inbreeding
temp <- sim_het %>%
  rename(sample = population) %>%
  group_by(sample) %>%
  summarise(het_mean = mean(het))

p2 <- phen_sum %>%
  filter(generation == 10) %>%
  mutate(sample = paste(nitrate_grown, selection, replicate, sep = "-")) %>%
  select(sample, f_mean, selection, nitrate_grown) %>%
  full_join(temp, by = "sample") %>%
  ggplot(aes(f_mean, het_mean, fill = selection)) +
  geom_point(shape = 21, size = 3) +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Inbreeding coefficient\n(Gen 10)",
       y = "Mean Heterozygosity\n(simulated)") +
  theme(legend.position = "none")

# plot the scan
p3 <- sim_het %>%
  # calculate mean and SD of simulated heterozygosity
  group_by(population) %>%
  summarise(sim_het_mean = mean(het),
            sim_het_sd = sd(het)) %>%
  # join with empirical data
  full_join(het, by = c("population" = "sample")) %>%
  filter(selection != "stabilising") %>%
  # scale as a z-score
  mutate(zscore = (het - sim_het_mean)/sim_het_sd) %>%
  # define a threshold as -4 SD below simulated means
  mutate(sig = ifelse(zscore <= -4, zscore, NA)) %>%
  mutate(selection = fct_relevel(selection, "random")) %>%
  ggplot(aes(pos/1e6, zscore)) +
  geom_line(aes(colour = selection)) +
  geom_point(aes(y = sig)) +
  geom_hline(yintercept = -4, linetype = "dashed", colour = "grey48") +
  geom_rect(data = centromeres, inherit.aes = FALSE,
            aes(xmin = start/1e6 - 0.5, xmax = end/1e6 + 0.5, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.5) +
  facet_grid(nitrate + rep ~ chrom, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = seq(0, 30, 10)) +
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "Mb", y = "Scaled heterozygosity") +
  theme(legend.position = "none")

pdf("./figures/FigS07.pdf", width = 7.5, height = 10)
(
  (
    (p1 | p2) + plot_layout(widths = c(2, 1))
  )
  / p3
) +
  plot_annotation(tag_levels = "A") +
  plot_layout(heights = c(1, 4), guides = "collect") &
  theme(legend.position = "top")
dev.off()
