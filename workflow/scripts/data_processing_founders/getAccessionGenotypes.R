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

if(!dir.exists(opt$outdir)){
  dir.create(opt$outdir)
  setwd(opt$outdir)
} else {
  setwd(opt$outdir)
}


#
# Read and tidy founder genotypes ----
#

# Link to server
server <- "http://mtweb.cs.ucl.ac.uk/mus/www/19genomes/variants.tables/"

# Download data for each chromosome
founder_genotypes <- map(1:5, function(i){
  # Name of the chromosome file
  chrom_file <- paste0("chr", i, ".alleles.txt")

  # Read the file - directly from web, takes a while
  chrom_data <- read_table2(file.path(server, chrom_file),
                            col_types = cols(.default = col_character(),
                                             pse.bp = col_integer(),
                                             bp = col_double(),
                                             nalleles = col_integer(),
                                             maf = col_integer()))

  # Tidy the file
  chrom_data <- chrom_data %>%
    # Rename some columns
    rename(chrom = chr, start = bp, major_allele_count = maf) %>%
    rename_all(list(~ str_replace(., "-.*", ""))) %>%
    # add end position - useful for finding overlaps later
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
    rename_at(vars(start), funs(deparse(i))) %>%
    rename(Ref = ref) %>%
    write_csv(paste0("accessions_chrom", i, ".dgrp"))

  # Return data in long format
  chrom_data %>%
    gather(accession, snp_allele, bur:zu)
})

# Bind the chromosome tables together
founder_genotypes <- bind_rows(founder_genotypes) %>%
  select(snp, chrom, start, end, major_allele_count, ref_allele_count, ref, accession, snp_allele)


#
# Get diagnostic SNPs ----
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


#
# Write tables ----
#

# List of founder SNP alleles including diagnostic accession if there is one
founder_genotypes %>%
  distinct(snp, chrom, start, end, ref, snp_allele, major_allele_count, ref_allele_count) %>%
  filter(ref != snp_allele) %>%
  left_join(diag_snps) %>%
  rename(alt = snp_allele) %>%
  write_csv("accessions_snps.csv")

# Write the genotypes themselves in tabular format
founder_genotypes %>%
  write_csv("accessions_genotypes.csv")
