##
# Script to process the phenotype data from selection experiment
# can be run non-interactively with RScript
##


#
# Setup ----
#
# Load packages
library(tidyverse)

# need `pedigree` package, not available through conda
if(!("pedigree" %in% installed.packages())){
  install.packages("pedigree")
}

# not loading the whole package due to conflicting functions
countGen <- pedigree::countGen
calcInbreeding <- pedigree::calcInbreeding


#
# Custom function ----
#
# This is a wrapper around `stats::wilcoxon.test()`, but outputs results in
# a tidy format with calculation of different effect sizes.
wilcoxonRankSum <- function(x, y, ...){
  # Remove missing values
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]

  # Sample sizes
  n_x <- length(x)
  n_y <- length(y)

  # Make the test and tidy it into a table
  tidy_test <- wilcox.test(x, y, ...) %>%
    broom::tidy()

  # Calculate effect sizes:
  # https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test#Effect_sizes

  # common language effect size (in both directions)
  x_less_y <- sum(as.vector(outer(x, y, "<")))/(n_x*n_y)
  y_less_x <- sum(as.vector(outer(x, y, ">")))/(n_x*n_y)
  x_equal_y <- sum(as.vector(outer(x, y, "==")))/(n_x*n_y)

  # Calculate rank-biserial correlation - this depends on the hypothesis
  # so using the output from the test
  r <- 1 - ((2*tidy_test$statistic)/(n_x*n_y))

  # Vargha and Delaney's A
  vda <- effsize::VD.A(x, y)$estimate

  # Output tidy data.frame
  tidy_test <- tidy_test %>%
    mutate(n_x = n_x,
           n_y = n_y,
           x_less_y = x_less_y,
           y_less_x = y_less_x,
           x_equal_y = x_equal_y,
           r = r,
           vda = vda)
}



#
# Read and filter data ----
#
# Read phenotype data
phen <- read_csv("./data/raw/phenotypes/compiled_phenotypes_data_clean.csv",
                 col_types = "iccccDiiinninccccccclnDiiiDnnncccllc")

# factorize date-related variables
phen <- phen %>%
  mutate(score_day = factor(score_day),
         score_month = factor(score_month),
         score_year = factor(score_year))

# Change selection names to be more informative
phen <- phen %>%
  mutate(selection = case_when(selection == "random" ~ "random",
                               selection == "most" ~ "directional",
                               selection == "average" ~ "stabilising",
                               TRUE ~ as.character(NA))) %>%
  mutate(selection = factor(selection, levels = c("random", "directional", "stabilising")))

# Some individuals of a family occur in different populations - might have been used in different crosses
# (except for generation 0, where this is expected)
# I left this unchanged, but something to be aware of.
phen %>%
  group_by(generation, family) %>%
  summarise(n_pop = length(unique(population)),
            pop = paste(unique(population), collapse = "; ")) %>%
  filter(n_pop > 1 & generation != 0) %>%
  as.data.frame()

# A related problem is that then the parent's population does not always match the progeny population
# These are mostly on generation 2 and 10
table(phen$generation, phen$population_discordant, useNA = "ifany")


# Keep only individuals that are not duplicated entries (only ~0.2%)
## also only keep individuals that were grown and selected in the same nitrate
## and remove the MAGIC line data from this analysis
phen_filter <- phen %>%
  filter(!duplicated_line | is.na(duplicated_line)) %>%
  filter(nitrate_grown == nitrate_selected | is.na(nitrate_selected)) %>%
  filter(generation != -1)



#
# Calculate inbreeding coefficients ----
#
# This confirms that the generations are being well interpreted by the pedigree package
temp <- phen_filter %>%
  select(line, mother, father) %>%
  as.data.frame() %>%
  countGen()

stopifnot(all(temp == phen_filter$generation))
rm(temp)


# Calculate inbreeding coefficients from pedigree
f <- phen_filter %>%
  select(line, mother, father) %>%
  as.data.frame() %>%
  calcInbreeding()

# Add to table
phen_filter$f <- f
rm(f)


#
# Summarise dataset ----
#
# 1 row per generation, nitrate, selection, replicate, population
# Calculate summary metrics and also selection differentials
phen_sum <- phen_filter %>%
  group_by(generation, nitrate_grown, selection, replicate, population) %>%
  rename(bolt = sow_to_bolt) %>%
  summarise_at(vars(total, height, bolt, f),
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
  group_by(population) %>%
  arrange(generation) %>%
  mutate(total_seldif = total_mean_sel - total_mean,
         height_seldif = height_mean_sel - height_mean,
         bolt_seldif = bolt_mean_sel - bolt_mean,
         total_cumsel = lag(cumsum(total_seldif)),
         height_cumsel = lag(cumsum(height_seldif)),
         bolt_cumsel = lag(cumsum(bolt_seldif)),
         total_sd_seldif = total_sd_sel - total_sd) %>%
  ungroup()




#
# Comparing control and selected populations ----
#
# Calculate difference between average of selected and random populations
phen_response <- phen_sum %>%
  group_by(generation, replicate, nitrate_grown) %>%
  summarise(cumsel_directional = total_cumsel[selection == "directional"],
            response_directional = total_mean[selection == "directional"] - total_mean[selection == "random"],
            cumsel_stabilising = total_cumsel[selection == "stabilising"],
            response_stabilising = total_mean[selection == "stabilising"] - total_mean[selection == "random"]) %>%
  ungroup() %>%
  gather(key, value, cumsel_directional:response_stabilising) %>%
  separate(key, c("metric", "comparison")) %>%
  spread(metric, value)


# Run non-parametric comparison and obtain effect sizes
## Comparing "directional" with "random"
directional_test <- phen_filter %>%
  mutate(id = 1:n()) %>%
  select(id, generation, replicate, nitrate_grown, total, selection) %>%
  drop_na() %>%
  spread(selection, total) %>%
  group_by(generation, nitrate_grown, replicate) %>%
  do(wilcoxonRankSum(.$directional, .$random, conf.int = TRUE)) %>%
  mutate(comparison = "directional") %>%
  ungroup()

## comparing "stabilising" with "random"
stabilising_test <- phen_filter %>%
  mutate(id = 1:n()) %>%
  select(id, generation, replicate, nitrate_grown, total, selection) %>%
  drop_na() %>%
  spread(selection, total) %>%
  group_by(generation, nitrate_grown, replicate) %>%
  do(wilcoxonRankSum(.$stabilising, .$random, conf.int = TRUE)) %>%
  mutate(comparison = "stabilising") %>%
  ungroup()

## bind the results in a single table
all_tests <- bind_rows(directional_test, stabilising_test) %>%
  ungroup() %>%
  mutate(p.adjust = p.adjust(p.value, method = "fdr"))

## Join with main response table
phen_response <- full_join(phen_response, all_tests,
          by = c("generation", "replicate", "nitrate_grown", "comparison"))

rm(directional_test, stabilising_test, all_tests)



#
# Comparing response in variance ----
#
# function to run test for equal variances
flignerTest <- function(x, y, ...){
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]

  # Output of fligner.test
  test_out <- fligner.test(list(x, y)) %>% broom::tidy()

  # Add difference of IQR
  test_out$iqr_dif <- IQR(y) - IQR(x)

  # Add difference of MAD
  test_out$mad_dif <- mad(y) - mad(x)

  # Add difference in SD
  test_out$sd_dif <- sd(y) - sd(x)

  # Add difference in variance
  test_out$var_dif <- var(y) - var(x)

  # Add difference in CV
  test_out$cv_dif <- sd(y)/mean(y) - sd(x)/mean(x)

  return(test_out)
}

directional_test <- phen_filter %>%
  mutate(id = 1:n()) %>%
  select(id, generation, replicate, nitrate_grown, total, selection) %>%
  drop_na() %>%
  spread(selection, total) %>%
  group_by(generation, nitrate_grown, replicate) %>%
  do(flignerTest(.$random, .$directional)) %>%
  mutate(comparison = "directional") %>%
  ungroup()

stabilising_test <- phen_filter %>%
  mutate(id = 1:n()) %>%
  select(id, generation, replicate, nitrate_grown, total, selection) %>%
  filter(complete.cases(.)) %>%
  spread(selection, total) %>%
  group_by(generation, nitrate_grown, replicate) %>%
  do(flignerTest(.$random, .$stabilising)) %>%
  mutate(comparison = "stabilising") %>%
  ungroup()

var_tests <- bind_rows(directional_test, stabilising_test) %>%
  ungroup() %>%
  mutate(p.adjust = p.adjust(p.value, method = "fdr"))

rm(directional_test, stabilising_test)


#
# Calculate plasticity in Gen 10 ----
#
gen10_plast <- phen %>%
  filter(generation == 10 & !duplicated_line & !population_discordant) %>%
  group_by(population, selection, nitrate_selected, replicate, family, mother, father) %>%
  rename(bolt = sow_to_bolt) %>%
  summarise_at(vars(total, bolt, height),
               list(mean_low = ~ mean(.[nitrate_grown == "LN"], na.rm = TRUE),
                    mean_high = ~ mean(.[nitrate_grown == "HN"], na.rm = TRUE),
                    mean_D = ~ mean(.[nitrate_grown == "HN"], na.rm = TRUE) - mean(.[nitrate_grown == "LN"], na.rm = TRUE),
                    var_low = ~ var(.[nitrate_grown == "LN"], na.rm = TRUE),
                    var_high = ~ var(.[nitrate_grown == "HN"], na.rm = TRUE),
                    n_low = ~ sum(!is.na(.[nitrate_grown == "LN"])),
                    n_high = ~ sum(!is.na(.[nitrate_grown == "HN"])))) %>%
  ungroup()


#
# Saving data ----
#
phen_filter %>%
  write_rds("./data/processed/phenotypes/phenotypes_individual.rds")

phen_sum %>%
  write_rds("./data/processed/phenotypes/phenotypes_summarised.rds")

phen_response %>%
  write_rds("./data/processed/phenotypes/phenotypes_response.rds")

var_tests %>%
  write_rds("./data/processed/phenotypes/variance_tests.rds")

gen10_plast %>%
  write_rds("./data/processed/phenotypes/generation10_plasticity.rds")





