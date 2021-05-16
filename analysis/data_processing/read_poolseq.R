##
# Script to read and process the poolseq data from selection experiment
##
library(tidyverse)

# read heterozygosity scans at 200Kb windows
pool200 <- read_csv("./data/processed/poolseq/haplotype_diversity_200kb.csv",
                    col_types = cols(
                      sample = col_character(),
                      nitrate = col_character(),
                      selection = col_character(),
                      rep = col_character(),
                      chrom = col_integer(),
                      pos = col_integer(),
                      start = col_integer(),
                      end = col_integer(),
                      hap_sum = col_double(),
                      hap_hom1 = col_double(),
                      hap_hom12 = col_double(),
                      hap_hom123 = col_double(),
                      hap_hom1234 = col_double(),
                      hap_hom2 = col_double(),
                      hap_shannon = col_double(),
                      hap_diversity = col_double(),
                      hap_nalleles = col_double(),
                      freqs_csv = col_character()
                    )) %>%
  mutate(selection = factor(selection, levels = c("random", "directional", "stabilising"))) %>%
  # add expected heterozygosity
  mutate(hap_het = 1 - hap_hom1)

# candidate sweeps - see "analysis/data_processing/identify_candidate_sweeps.R"
candidate_sweeps <- read_csv("./data/processed/poolseq/candidate_sweep_intervals.csv",
                               col_types = cols())


# Read centromere region as annotated in TAIR9 release
# ftp://ftp.arabidopsis.org/home/tair/Genes/TAIR9_genome_release/TAIR9_gff3/Assembly_GFF
# see https://www.biostars.org/p/18782/
centromeres <- read_csv("./data/external/tair9_centromeres.csv",
                        col_types = cols())


# Define thresholds for peak identification of each heterozygosity statistic
threshold <- pool200 %>%
  filter(selection == "directional") %>%
  group_by(sample) %>%
  summarise(across(matches("hap_hom1"),
                   list(q1 = ~ quantile(1 - ., 0.01)))) %>%
  ungroup() %>%
  summarise(across(matches("_q1"), median))


