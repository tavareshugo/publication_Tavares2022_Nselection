library(tidyverse)

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

# read processed phenotype data to fetch generation 10 individuals in each population
gen10_idx <- read_rds("data/processed/phenotypes/phenotypes_individual.rds") %>%
  select(generation, line, nitrate_grown, selection, replicate,
         line, mother, father) %>%
  filter(generation == 10) %>%
  mutate(population = paste(nitrate_grown, selection, replicate, sep = "-")) %>%
  with(split(line, population))


# Simulate Genotypes --------------------------------------------------------------

# simulate genotype frequency change
simGenoFreqChange <- function(start_freq = 1/19){
  message("Starting frequency = ", start_freq, "\nseed = ", round(1596467939*start_freq))
  set.seed(round(1596467939*start_freq)) # for reproducibility
  
  # initiate a matrix of genotypes
  n_sim_loci <- 1000 # number of loci to simulate
  genotypes <- matrix(nrow = nrow(ped), ncol = n_sim_loci*2)
  rownames(genotypes) <- ped$line
  colnames(genotypes) <- paste0("locus", rep(1:n_sim_loci, each = 2), ".", c(1, 2))

  # function to simulate starting genotypes (for the parents)
  simStartGeno <- function(nloci = n_sim_loci, start_freq){
    out <- rbinom(nloci, size = 1, prob = start_freq)
    out <- rep(out, each = 2) # homozygous
    return(out)
  }

  for(i in 1:nrow(ped)){
    # get IDs of the current individual and its parents
    id <- ped$line[i]
    mother <- ped$mother[i]
    father <- ped$father[i]

    # cat(i, id, mother, father, "\n") # print message for tracking progress

    # sometimes parents are repeated so we should only sample them once
    if(any(is.na(genotypes[id, ]))){
      if(any(!is.na(genotypes[id, ]))) stop("Some NA's!?")

      # simulate genotype if it's the founder generation
      if(ped$generation[i] == -1){
        genotypes[id, ] <- simStartGeno(nloci = n_sim_loci, start_freq = start_freq)

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

  # clear workspace
  rm(i, id, mother, father, maternal_alleles, paternal_alleles)

  # Calculate change in frequency ------------------------

  # calculate allele frequencies in generation 10
  sim_freqs <- lapply(gen10_idx, function(ids){

    # sample from the genotypes to account for PoolSeq sampling step
    # assuming depth of coverage 100x
    sampled_ids <- sample(ids, 100, replace = TRUE)

    out <- colMeans(genotypes[sampled_ids, ])
    names(out) <- colnames(genotypes)

    return(out)
  })

  # tidy things into a data.frame
  sim_freqs <- map_dfr(sim_freqs, function(d){
    d %>%
      enframe(name = "locus", value = "freq")
  }, .id = "population")

  # summarise the two copies of each locus as one
  sim_freqs <- sim_freqs %>%
    mutate(locus = str_remove(locus, "locus")) %>%
    separate(locus, c("locus", "copy"), convert = TRUE) %>%
    group_by(population, locus) %>%
    summarise(freq = sum(freq)/2) %>%
    ungroup()
    
  return(sim_freqs)

}

# simulate frequencies for different starting frequencies
sim_freqs <- map((1:9)/19, simGenoFreqChange)
names(sim_freqs) <- 1:9
sim_freqs <- bind_rows(sim_freqs, .id = "start_mac")

write_csv(sim_freqs, 
          "./data/processed/snp_analysis/simulated_frequency_changes.csv")


# Estimate Thresholds ------

# arcsin transform
sqrtasin <- function(x){
  2*asin(sqrt(x))
}

sim_freqs %>%
  mutate(start_mac = as.numeric(start_mac)) %>%
  # calculate squared arcsin frequency difference - like Castro et al.
  # and absolute allele frequency change
  mutate(z2 = (sqrtasin(freq) - sqrtasin(start_mac/19))^2 / pi^2,
         afc = abs(freq - start_mac/19)) %>% 
  # estimate summary stats
  group_by(start_mac, population) %>%
  summarise(z2_q99 = quantile(z2, 0.99),
            afc_q99 = quantile(afc, 0.99),
            z2_mean = mean(z2),
            z2_sd = sd(z2),
            afc_mean = mean(afc),
            afc_sd = sd(afc)) %>%
  ungroup() %>%
  write_csv("./data/processed/snp_analysis/simulated_thresholds.csv")
