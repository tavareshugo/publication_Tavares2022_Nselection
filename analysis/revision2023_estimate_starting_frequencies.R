library(tidyverse)

# Read pedigree -----------------------------------------------------------

# get MAGIC/HSRIL ids
ids <- read_csv("./data/external/founder_genotypes/magic_mosaics.csv") %>%
  distinct(hsril, magic)

# get experimental pedigree
ped <- read_csv("./data/raw/phenotypes/compiled_phenotypes_data_clean.csv",
                col_types = "iccccDiiinninccccccclnDiiiDnnncccllc")

# remove few individuals without complete pedigree information
ped <- ped %>%
  mutate(population = paste(nitrate_selected, selection, replicate, sep = "-")) %>%
  mutate(population = population %>%
                        str_replace("most", "directional") %>%
                        str_replace("average", "stabilising")) %>%
  distinct(population, generation, line, mother, father)

# there are some founders with missing phenotype data, so we add them separately
miss <- ped %>%
  filter(!(father %in% line & mother %in% line) & generation != -1) %>%
  select(population, mother, father) %>%
  pivot_longer(c(mother, father), values_to = "line") %>%
  distinct(line, population) %>%
  mutate(generation = -1L, mother = NA, father = NA) %>%
  select(generation, line, population, mother, father)

ped <- bind_rows(miss, ped) %>%
  distinct(generation, line, population, mother, father) %>%
  arrange(population, generation, line, mother, father)

# get the founder IDs
founders <- ped %>%
  filter(generation == -1) %>%
  select(line, population) %>%
  group_by(population) %>%
  group_nest() %>%
  mutate(data = map(data, ~ .x$line))

rm(miss, ped)


# Read genotypes --------

# readr::read_csv doesn't work so well
magic_snps <- data.table::fread("./data/external/founder_genotypes/imputed_magic_snp.csv")
magic_snps <- magic_snps %>%
  as.data.frame() %>%
  column_to_rownames("magic") %>%
  as.matrix()

# there are too many MAGIC lines with missing genotype data
# they have not been sequenced
founders %>%
  mutate(miss = map_dbl(data, ~ sum(!.x %in% ids$hsril)),
         total = map_dbl(data, length)) %>%
  mutate(fraction_miss = miss/total)

# therefore we chose to approximate the starting frequencies
# by assuming each accession starts at ~1/19 frequency


# if we were being perfectionist we should try to impute all those genotypes
# as there will be a deviation from the 1/19 expectation when we have 40 founders
sim <- lapply(1:1000, function(i){
  tibble(allele = sample(letters[1:19], size = 40*2, replace = TRUE)) |> 
    count(allele) |> 
    mutate(rep = i)
}) |> 
  bind_rows()
sim |> 
  complete(rep, allele) |> 
  mutate(n = replace_na(n, 0)) |> 
  group_by(rep) |> 
  mutate(freq = n/sum(n)) |> 
  ungroup() |> 
  ggplot(aes(freq, group = allele)) +
  geom_freqpoly(binwidth = 1/80) + 
  geom_vline(xintercept = 1/19) +
  scale_x_continuous(limits = c(0, 0.3))
