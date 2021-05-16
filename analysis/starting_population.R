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


# Get MAGIC genotypes ----------------------------------------------------

# Load the data to the environment
data("kover2009")

# impute most likely genotypes
magics <- viterbi(kover2009, kover2009$gmap, error_prob = 0.01)
magics <- as.matrix(magics)
magics <- tolower(magics)

# get public names
tf <- tempfile(fileext = ".xls")
download.file("https://www.genetics.org/lookup/suppl/doi:10.1534/genetics.114.170746/-/DC1/genetics.114.170746-3.xls", tf)
magic_names <- readxl::read_xlsx(tf,
                                 sheet = "seed area", range = "H1:I824")

rownames(magics) <- magic_names$`our name`[match(rownames(magics), magic_names$`Public name`)]

# retain only those occurring in our data
magics <- magics[(rownames(magics) %in% ped$line), ]

# clean workspace
unlink(tf); rm(tf, magic_names, kover2009)


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

freqs <- apply(magics, 2, function(i) enframe(table(i)))
freqs <- bind_rows(freqs, .id = "locus")
freqs %>%
  group_by(locus) %>%
  mutate(freq = value/sum(value)) %>%
  ungroup() %>%
  ggplot(aes(freq, name)) +
  ggridges::geom_density_ridges(fill = "dodgerblue", alpha = 0.5) +
  geom_vline(xintercept = 1/19, linetype = "dashed") +
  annotate(geom = "label", x = 1/19, y = -1, label = "1/19", vjust = 0) +
  labs(x = "Frequency", y = "Accession", title = "Allele frequency in MAGIC lines")
