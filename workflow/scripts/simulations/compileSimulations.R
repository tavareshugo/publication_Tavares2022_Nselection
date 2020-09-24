library(tidyverse)

# custom read function
read_sim <- function(prefix, n_selected_loci, selected_effect,
                     n_adv_alleles,
                     n = 5){
  # list input files
  infiles <- list.files("data/intermediate/simulations/",
             pattern = paste0(prefix, n_selected_loci, "-",
                              selected_effect, "-", n_adv_alleles),
             full.names = TRUE)
  names(infiles) <- infiles %>% basename() %>%
    str_remove(prefix) %>% str_remove(".csv")

  if(length(infiles) == 0) stop("no files found: ", paste0(prefix, n_selected_loci, "-", selected_effect, "-", n_adv_alleles))

  if(n < 1 | n > length(infiles)) n <- length(infiles)

  # set column types for each type of file (for faster reading)
  if(prefix == "het_"){
    coltypes = "ciclnnnnnnnn"
  } else if (prefix == "pedigree_"){
    coltypes = "cicccnnnn"
  } else if (prefix == "freq_"){
    coltypes = "cicnl"
  } else {
    stop("prefix has to be one of: het_, freq_, pedigree_")
  }

  # read data into single data frame
  x <- map_dfr(infiles[1:n], read_csv, col_types = coltypes, .id = "id")
  x
}

# read parameter file list
params <- read_csv("data/intermediate/simulations/parameter_list.csv")


# Heterozygosity ------------------------------------------------------
message("Processing heterozygosity files...")

# empty lists to store outputs
gen10_het <- vector("list", length = nrow(params))
gen10_het_summary <- vector("list", length = nrow(params))
gen10_het_selected_loci <- vector("list", length = nrow(params))

# loop through parameter file
for (i in 1:nrow(params)){
  n_selected_loci <- params$n_selected_loci[i]
  selected_effect <- params$selected_effect[i]
  n_adv_alleles <- params$n_adv_alleles[i]
  output_suffix <- paste(n_selected_loci, selected_effect, n_adv_alleles, sep = "-")

  message("Processing ", output_suffix)

  het <- read_sim("het_",
                  n_selected_loci = n_selected_loci,
                  selected_effect = selected_effect,
                  n_adv_alleles = n_adv_alleles,
                  n = Inf)

  # save to list
  gen10_het[[i]] <- het %>%
    filter(gen == 10) %>%
    separate(locus, into = c("chrom", "pos"), sep = "-")

  # get heterozygosity for selected loci only
  gen10_het_selected_loci[[i]] <- gen10_het[[i]] %>%
    filter(selected)

  # summarise across replicates
  gen10_het_summary[[i]] <- gen10_het[[i]] %>%
    group_by(selection, gen, locus, selected) %>%
    summarise(across(starts_with("het"),
                     list(mean = mean,
                          median = median,
                          min = min,
                          max = max,
                          q10 = ~ quantile(., probs = 0.1),
                          q90 = ~ quantile(., probs = 0.9)))) %>%
    ungroup() %>%
    separate(locus, into = c("chrom", "pos"), sep = "-")
}

# bind and write results
gen10_het %>%
  bind_rows() %>%
  write_csv("data/processed/simulations/het_all_loci.csv")

gen10_het_summary %>%
  bind_rows() %>%
  write_csv("data/processed/simulations/het_summarised.csv")

gen10_het_selected_loci %>%
  bind_rows() %>%
  write_csv("data/processed/simulations/het_selected_loci.csv")

# gen10_het_selected_loci[[1]] %>%
#   pivot_longer(starts_with("het")) %>%
#   ggplot(aes(locus, value, fill = selection)) +
#   ggbeeswarm::geom_quasirandom(aes(colour = selection),
#                                dodge.width = 0.5) +
#   facet_wrap(~ name)
#
# gen10_het[[1]] %>%
#   separate(locus, into = c("chrom", "pos"), sep = "-") %>%
#   ggplot(aes(as.integer(pos), het)) +
#   geom_line(aes(group = interaction(selection, id)), alpha = 0.1) +
#   theme_void() +
#   facet_grid(selection ~ chrom, scale = "free_x", space = "free_x")
#
# gen10_het_summary[[1]] %>%
#   ggplot(aes(as.integer(pos), het_mean)) +
#   geom_ribbon(aes(ymin = het_q10, ymax = het_q90), alpha = 0.1) +
#   geom_line() +
#   theme_void() +
#   scale_colour_viridis_d() +
#   facet_grid(selection ~ chrom, scale = "free_x", space = "free_x")


# Pedigree ----------------------------------------------------------------
message("Processing train/pedigree files...")

# function to calculate inbreeding
inbreeding <- function(ind, mother, father){
  d <- data.frame(ind, mother, father)
  pedigree::calcInbreeding(d)
}

# empty lists to store outputs
ped_summary <- vector("list", length = nrow(params))
trait_response <- vector("list", length = nrow(params))
trait_h2 <- vector("list", length = nrow(params))

# loop through parameter file
for (i in 1:nrow(params)){
  n_selected_loci <- params$n_selected_loci[i]
  selected_effect <- params$selected_effect[i]
  n_adv_alleles <- params$n_adv_alleles[i]
  output_suffix <- paste(n_selected_loci, selected_effect, n_adv_alleles, sep = "-")

  message("Processing ", output_suffix)

  ped <- read_sim("pedigree_",
                  n_selected_loci = n_selected_loci,
                  selected_effect = selected_effect,
                  n_adv_alleles = n_adv_alleles,
                  n = Inf)

  # save to list
  ped_summary[[i]] <- ped %>%
    group_by(id, selection) %>%
    mutate(f = inbreeding(ind_id, mother_id, father_id)) %>%
    group_by(id, selection, gen) %>%
    summarise(across(c("trait", "f"),
                     list(mean = mean,
                          median = median,
                          min = min,
                          max = max,
                          sd = sd,
                          q10 = ~ quantile(., probs = 0.1),
                          q90 = ~ quantile(., probs = 0.9),
                          mean_sel = ~ mean(.[fitness == 1]),
                          sd_sel = ~ sd(.[fitness == 1])))) %>%
    group_by(id, selection) %>%
    arrange(gen) %>%
    mutate(seldif = trait_mean_sel - trait_mean,
           cumsel = lag(cumsum(seldif)),
           selintensity = seldif/trait_sd) %>%
    ungroup()

  trait_response[[i]] <- ped_summary[[i]] %>%
    group_by(id, gen) %>%
    summarise(cumsel_directional = cumsel[selection == "directional"],
              response_directional = trait_mean[selection == "directional"] - trait_mean[selection == "random"],
              cumsel_stabilising = cumsel[selection == "stabilising"],
              response_stabilising = trait_mean[selection == "stabilising"] - trait_mean[selection == "random"]) %>%
    ungroup() %>%
    gather(key, value, cumsel_directional:response_stabilising) %>%
    separate(key, c("metric", "comparison")) %>%
    spread(metric, value)

  trait_h2[[i]] <- trait_response[[i]] %>%
    filter(comparison == "directional") %>%
    group_by(id) %>%
    nest() %>%
    mutate(h2 = map(data, function(i){
      lm(response ~ 0 + cumsel, data = i) %>%
        broom::tidy()
  })) %>%
    select(-data) %>%
    unnest(h2)

}

# bind and export
ped_summary %>%
  bind_rows() %>%
  write_csv("data/processed/simulations/trait_summarised.csv")

trait_response %>%
  bind_rows() %>%
  write_csv("data/processed/simulations/trait_responses.csv")

trait_h2 %>%
  bind_rows() %>%
  write_csv("data/processed/simulations/trait_heritabilities.csv")


# ped_summary[[i]] %>%
#   ggplot(aes(gen, f_mean)) +
#   geom_ribbon(stat = "summary", fun.data = "mean_se",
#               fun.args = list(mult = 1.96),
#               aes(fill = selection), alpha = 0.5)
#
# trait_response[[i]] %>%
#   filter(comparison == "directional") %>%
#   ggplot(aes(cumsel, response)) +
#   geom_line(aes(group = id), alpha = 0.5, colour = "grey") +
#   geom_smooth(se = FALSE, colour = "black")
#
# trait_h2[[i]] %>%
#   ggplot(aes(estimate)) +
#   geom_density()
