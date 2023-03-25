##
# Script to read and process the snape data from selection experiment
##
library(tidyverse)
library(windowscanr) # https://github.com/tavareshugo/WindowScanR

# arcsin transform
sqrtasin <- function(x){
  2*asin(sqrt(x))
}


# all population files
snape <- tibble(file = c("hn_average_a_reseq.snape.csv",
           "hn_average_b_reseq.snape.csv",
           "hn_average_c.snape.csv",
           "hn_most_a.snape.csv",
           "hn_most_b_reseq.snape.csv",
           "hn_most_c.snape.csv",
           "hn_random_a.snape.csv",
           "hn_random_b.snape.csv",
           "hn_random_c.snape.csv",
           "ln_average_a.snape.csv",
           "ln_average_b_reseq.snape.csv",
           "ln_average_c.snape.csv",
           "ln_most_a_reseq.snape.csv",
           "ln_most_b.snape.csv",
           "ln_most_c.snape.csv",
           "ln_random_a_reseq.snape.csv",
           "ln_random_b.snape.csv",
           "ln_random_c.snape.csv"))

# table
snape <- snape %>%
	separate(file, c("nitrate", "selection", "rep"), sep = "_", remove = FALSE) %>%
  mutate(rep = str_remove(rep, ".snape.csv")) %>%
  mutate(across(c(nitrate, rep), toupper))

snape <- snape %>%
  mutate(data = map(file, ~ read_csv(paste0("./data/revisions/snape/filtered/", .x), 
                                     col_types = "iiciiiicddd")))

# unnest table and tidy
snape <- snape %>%
  select(-file) %>%
  unnest(data)  %>%
  # tidy alleles
  rename(ref = ref_allele)  %>%
  mutate(ref = toupper(ref), alleles = toupper(alleles)) %>%
  # rename selection variable
  mutate(selection = case_when(selection == "most" ~ "directional",
                               selection == "average" ~ "stabilising",
                               selection == "random" ~ "random",
                               TRUE ~ NA_character_))


# Calculate AFC -----

# parent snps
snps <- read_csv("./data/external/founder_genotypes/accessions_snps.csv")
snps <- snps |> 
  select(chrom, pos = start,
         start_ref_allele_count = ref_allele_count, 
         ref, alt)

# join to snape frequencies
snape <- inner_join(snape, snps, by = c("chrom", "pos", "ref"))

# calculate AFC
snape <- snape %>%
  # express frequency relative to the starting minor allele
  mutate(maf_end = ifelse(start_ref_allele_count > 9, 
                          alt_freq, 
                          1 - alt_freq),
         mac_start = ifelse(start_ref_allele_count > 9,
                            19 - start_ref_allele_count,
                            start_ref_allele_count)) %>%
  mutate(maf_start = mac_start/19) %>%
  mutate(z2 = (sqrtasin(maf_end) - sqrtasin(maf_start))^2 / pi^2,
         afc = abs(maf_end - maf_start))

rm(snps)

# Add thresholds -----

# frequency change thresholds from simulations
thr <- read_csv("./data/processed/snp_analysis/simulated_thresholds.csv")

# add thresholds
snape <- snape %>%
  mutate(population = paste(nitrate, selection, rep, sep = "-")) %>%
  left_join(thr, by = c("population", "mac_start" = "start_mac"))

rm(thr)


# Summarise across windows -----


snape200 <- snape %>%
  mutate(z2_outlier = z2 > z2_q99, 
         z2_zscore = (z2 - z2_mean)/z2_sd,
         afc_outlier = afc > afc_q99,
         afc_zscore = (afc - afc_mean)/afc_sd) %>%
  winScan(groups = c("nitrate", "selection", "rep", "chrom"),
          position = "pos",
          values = c("z2_outlier", "afc_outlier", "afc", "z2", "afc_zscore", "z2_zscore"),
          win_size = 200e3, 
          win_step = 100e3, 
          funs = c("mean"), 
          cores = 4) %>%
  as_tibble()

write_csv(snape200, "./data/processed/poolseq/snp_freq_200kb.csv")
