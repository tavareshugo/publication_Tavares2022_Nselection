library(atMAGIC) # this loads qtl2 as well
library(qtl2helper)
library(tidyverse)
library(vroom)
theme_set(theme_classic())


# Get MAGIC founders ------------------------------------------------------

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


# Get window positions -----

# this will be used as a map for inferring accession genotypes
# to match window positions with the poolseq data
map <- read_csv("./data/processed/poolseq/haplotype_diversity_200kb.csv") %>%
  distinct(chrom, pos) %>%
  with(split(pos, chrom))

# the genetic map positions in the MAGIC lines are just a linear transformation
# of the physical map - we apply the same transformation here
# lm(unlist(kover2009$pmap) ~ unlist(kover2009$gmap))
map <- lapply(map, function(i) i/239884)
map <- insert_pseudomarkers(kover2009$gmap, pseudomarker_map = map)


# Get MAGIC genotypes ----------------------------------------------------

# Load the data to the environment
data("kover2009")

# impute most likely genotypes
#magics <- viterbi(kover2009, kover2009$gmap, error_prob = 0.01)
magics <- viterbi(kover2009, map, error_prob = 0.01)
magics <- as.matrix(magics)
magics <- tolower(magics)

# get public names
magic_names <- readxl::read_xlsx("./data/external/founder_genotypes/genetics.114.170746-3.xls",
                                 sheet = "seed area", range = "H1:I824")

rownames(magics) <- magic_names$`our name`[match(rownames(magics), magic_names$`Public name`)]

# retain only those occurring in our data
magics <- magics[(rownames(magics) %in% ped$line), ]


# Get accession genotypes -------------------------------------------------

acc <- vroom("data/external/founder_genotypes/accessions_genotypes.csv")

acc %>%
  # get minor alele counts
  distinct(snp, major_allele_count) %>%
  mutate(mac = 19 - major_allele_count) %>%
  # calculate their frequency
  count(mac) %>%
  mutate(frequency = n/sum(n)) %>%
  # plot
  rename(`Minor allele count` = mac) %>%
  ggplot(aes(`Minor allele count`, frequency)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  scale_x_continuous(breaks = 1:19)


# Calculate accession frequency -------------------------------------------

# counts each accession at each marker
freqs <- apply(magics, 2, function(i) enframe(table(i)))
freqs <- bind_rows(freqs, .id = "locus") %>%
  mutate(value = as.numeric(value))

# join with map position
freqs <- lapply(map, enframe, name = "locus", value = "pos") %>%
  bind_rows(.id = "chrom") %>%
  # get only the window markers
  filter(str_detect(locus, "^c")) %>%
  # transform back to physical distance
  mutate(pos = as.integer(pos*239884)) %>%
  inner_join(freqs, by = "locus") %>%
  select(-locus)

# calculate frequency
freqs <- freqs %>%
  group_by(chrom, pos) %>%
  mutate(freq = value/sum(value)) %>%
  ungroup()

freqs_sum <- freqs %>% 
  group_by(chrom, pos) %>%
  summarise(hap_het = 1-sum(freq^2),
            hap_hom1 = sum(freq^2),
            hap_hom12 = hap_hom1 + 2*sort(freq, decreasing = TRUE)[1]*sort(freq, decreasing = TRUE)[2],
            hap_hom2 = hap_hom1 - sort(freq, decreasing = TRUE)[1]^2,
            hap_shannon = -sum(freq[freq > 0] * log(freq[freq > 0])),
            hap_diversity = exp(hap_shannon),
            acc_freq_csv = paste(name, freq, sep = ",", collapse = "\n"))

write_csv(freqs_sum, "./data/external/founder_genotypes/startpop_window200kb.csv")


# visualise
freqs %>%
  ggplot(aes(freq, name)) +
  ggridges::geom_density_ridges(fill = "dodgerblue", alpha = 0.5) +
  geom_vline(xintercept = 1/19, linetype = "dashed") +
  annotate(geom = "label", x = 1/19, y = -1, label = "1/19", vjust = 0) +
  labs(x = "Frequency", y = "Accession", title = "Allele frequency in MAGIC lines")
