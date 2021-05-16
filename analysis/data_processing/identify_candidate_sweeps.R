#
# Identify candidate sweeps
#

library(tidyverse)
library(valr)

# read custom functions
source("./analysis/functions/findPeaks.R")



#### Read data ####

# heterozygosity scans at 200Kb windows
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


#### Define sweep identifier function ####

# custom function to find a sweep based on a particular homozygosity statistic
find_sweep <- function(statistic = "hap_hom1", x){
  require(valr)

  # make column with heterozygosity of interest
  x$het <- 1 - x[[statistic]]

  # Define threshold
  threshold <- x %>%
    filter(selection == "directional") %>%
    group_by(sample) %>%
    summarise(q1 = quantile(het, 0.01)) %>%
    summarise(threshold = median(q1)) %>%
    pull(threshold)
  print(threshold)

  # find sweeps
  candidate_sweeps <- x %>%
    #filter(selection != "stabilising") %>%
    group_by(sample) %>%
    mutate(peak = findPeaks(het, threshold, below = TRUE)) %>%
    drop_na(peak) %>%
    group_by(sample, nitrate, rep, selection, chrom, peak,
             .drop = TRUE) %>%
    summarise(start = min(start),
              end = max(end),
              pos = pos[which(het == min(het))],
              het = het[which(het == min(het))]) %>%
    ungroup()

  # merge intervals that are within 400Kb of each other
  candidate_sweeps <- candidate_sweeps %>%
    group_by(sample, nitrate, rep, selection, chrom) %>%
    bed_cluster(max_dist = 400e3) %>%
    group_by(sample, nitrate, rep, selection, chrom, .id) %>%
    summarise(start = min(start), end = max(end), pos = pos[which(het == min(het))]) %>%
    group_by(selection) %>%
    arrange(chrom, start) %>%
    mutate(peak_id = 1:n()) %>%
    ungroup() %>%
    select(-.id) %>%
    mutate(statistic = paste0("H", str_remove(statistic, "hap_hom")))

  return(candidate_sweeps)

}


#### identify sweeps ####

# run function for each statistic
candidate_sweeps <- map_df(c("hap_hom1", "hap_hom12", "hap_hom123", "hap_hom1234"),
                           find_sweep, x = pool200)

# write result
candidate_sweeps %>%
  write_csv("./data/processed/poolseq/candidate_sweep_intervals.csv")

#### deprecated ####

# #### identify peaks ####
#
# # Define thresholds for peak identification of each heterozygosity statistic
# threshold <- pool200 %>%
#   filter(selection == "directional") %>%
#   group_by(sample) %>%
#   summarise_at(vars(matches("hap_hom1")),
#                list(q1 = ~ quantile(1 - ., 0.01))) %>%
#   summarise_at(vars(matches("_q1")), median)
#
# # Find dips in heterozygosity
# candidate_sweeps <- pool200 %>%
#   filter(selection != "stabilising") %>%
#   group_by(sample) %>%
#   mutate(peak_h1 = findPeaks(1 - hap_hom1, threshold$hap_hom1_q1, below = TRUE),
#          peak_h12 = findPeaks(1 - hap_hom12, threshold$hap_hom12_q1, below = TRUE),
#          peak_h123 = findPeaks(1 - hap_hom123, threshold$hap_hom123_q1, below = TRUE),
#          peak_h1234 = findPeaks(1 - hap_hom1234, threshold$hap_hom1234_q1, below = TRUE)) %>%
#   gather("statistic", "peak", matches("peak_")) %>%
#   drop_na(peak) %>%
#   group_by(sample, nitrate, rep, selection, chrom, statistic, peak,
#            .drop = TRUE) %>%
#   summarise(start = min(start),
#             end = max(end)) %>%
#   ungroup() %>%
#   mutate(statistic = toupper(str_remove(statistic, "peak_")))
#
#
#
# # merge intervals that are within 400Kb of each other
# candidate_sweeps <- candidate_sweeps %>%
#   group_by(sample, nitrate, rep, selection, chrom, statistic) %>%
#   bed_cluster(max_dist = 400e3) %>%
#   group_by(sample, nitrate, rep, selection, chrom, statistic, .id) %>%
#   summarise(start = min(start), end = max(end)) %>%
#   group_by(statistic) %>%
#   arrange(desc(selection), chrom, start) %>%
#   mutate(peak_id = 1:n()) %>%
#   ungroup() %>%
#   select(-.id)
#
#
# # write results
# candidate_sweeps %>%
#   write_csv("./data/processed/candidate_sweep_intervals.csv")
