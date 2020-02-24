##################################################
# Prepare data with start population frequencies #
##################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(optparse)
})

option_list = list(
  make_option(c("--outdir"),
              action = "store",
              default = NA,
              type = 'character',
              help = "The directory to store output files.")
  )

opt <- parse_args(OptionParser(option_list=option_list))

opt$outdir <- "data/external/snps/"

if(!dir.exists(opt$outdir)){
  dir.create(opt$outdir)
  setwd(opt$outdir)
} else {
  setwd(opt$outdir)
}



#
# Read and tidy founder genotypes ----
#
## Download founder genotype data
## Read it and filter/tidy it

# Get data for each chromosome
founder_genotypes <- map(1:5, function(i){

  # server file
  remote_file <- paste0("http://mtweb.cs.ucl.ac.uk/mus/www/19genomes/variants.tables/chr", i, ".alleles.txt")

  # Read the file
  message("Downloading: ", remote_file)
  chrom_data <- read_table2(remote_file,
                            col_types = cols(.default = col_character(),
                                             pse.bp = col_integer(),
                                             bp = col_double(),
                                             nalleles = col_integer(),
                                             maf = col_integer()))

  # Tidy the file
  chrom_data <- chrom_data %>%
    # Rename some columns
    rename(chrom = chr, start = bp, major_allele_count = maf) %>%
    rename_all(funs(str_replace(., "-.*", ""))) %>%
    # add end position
    mutate(start = as.integer(start), end = start,
           ref = col, snp = paste(chrom, start, sep = "_")) %>%
    # filter all genotypes to retain only simple SNPs
    filter_at(vars(bur:zu), all_vars(. %in% c("A", "T", "C", "G"))) %>%
    # retain only bi-allelic sites
    filter(nalleles == 2 & start %% 1 == 0) %>%
    # get reference allele count (easier to detect diagnostic alleles)
    unite("alleles", bur:zu, sep = " ", remove = FALSE) %>%
    mutate(ref_allele_count = str_count(alleles, ref)) %>%
    select(-alleles)

  # Write dgrp formatted tables (for HARP software)
  chrom_data %>%
    select(start, ref, bur:zu) %>%
    rename_at(vars(start), funs(deparse(as.numeric(i)))) %>%
    rename(Ref = ref) %>%
    write_csv(paste0("founders_chrom", i, ".dgrp"))

  # Return data in long format
  chrom_data %>%
    gather(accession, snp_allele, bur:zu)
})

# Bind the chromosome tables together
founder_genotypes <- bind_rows(founder_genotypes) %>%
  select(snp, chrom, start, end, major_allele_count, ref_allele_count, ref, accession, snp_allele)


#
# Export SNPs ----
#
# Classify SNPs as diagnostic
founder_genotypes <- founder_genotypes %>%
  mutate(diag_accession = case_when(major_allele_count == 18 & ref_allele_count == 1 & ref == snp_allele ~ accession,
                                    major_allele_count == 18 & ref_allele_count == 18 & ref != snp_allele ~ accession,
                                    TRUE ~ as.character(NA)))

# Make table of diagnostic snps
diag_snps <- founder_genotypes %>%
  filter(!is.na(diag_accession)) %>%
  select(snp, ref, snp_allele, diag_accession)

# Write a list of founder SNP alleles including diagnostic accession if there is one
founder_genotypes %>%
  distinct(snp, chrom, start, end, ref, snp_allele) %>%
  filter(ref != snp_allele) %>%
  left_join(diag_snps) %>%
  rename(alt = snp_allele) %>%
  write_tsv("founders_snps.tsv")

# Write the genotypes themselves in tabular format
founder_genotypes %>%
  write_tsv("founders_genotypes.tsv")

