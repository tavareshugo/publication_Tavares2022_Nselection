############################
# Compile pool-seq metrics #
############################

suppressPackageStartupMessages({
  library(tidyverse)
  library(optparse)
})

option_list = list(
  make_option("--indir",
              action = "store",
              default = NA,
              type = 'character',
              help = "The directory with frequency files."),
  make_option("--outdir",
              action = "store",
              default = NA,
              type = 'character',
              help = "The directory to store output files."),
  make_option("--window",
              action = "store",
              default = NA,
              type = 'character',
              help = "Window size.")
)

opt <- parse_args(OptionParser(option_list=option_list))

message("outdir: ", opt$outdir)
message("indir: ", opt$indir)
message("window: ", opt$window)
# opt$outdir <- "data/processed/poolseq/"
# opt$indir <- "data/intermediate/harp_freq/"
# opt$window <- "200"

if(!dir.exists(opt$outdir)){
  dir.create(opt$outdir, recursive = TRUE)
}


#
# Custom functions ----
#
#' Read HARP output frequency files
#' @param pattern the pattern to match file names to read in
#' @param dir the directory to read files from
readHarp <- function(pattern, dir = "."){
  # List files based on pattern
  harp_files <- list.files(dir, pattern, full.names = TRUE)

  # Read all files into a named list
  harp_freqs <- map(harp_files, read.table,
                    header = FALSE, stringsAsFactors = FALSE,
                    col.names = c("chrom", "start", "end", "bur", "can", "col", "ct", "edi", "hi", "kn", "ler", "mt", "no", "oy", "po", "rsch", "sf", "tsu", "wil", "ws", "wu", "zu"))
  names(harp_freqs) <- basename(harp_files) %>% gsub("\\..*", "", .)

  # Bind the files to a single data.frame and convert to long format
  # note there's a warning from separate, which is fine (related to resequenced samples)
  harp_freqs <- harp_freqs %>%
    bind_rows(.id = "sample") %>%
    mutate(sample = str_remove(sample, "_chrom.*kb$")) %>%
    gather(acc, freq, bur:zu) %>%
    mutate(pos = (start + end)/2) %>%
    separate(sample, c("nitrate", "selection", "rep"),
             sep = "-", remove = FALSE)

  return(harp_freqs)
}


#' Estimate expected homozigosity
#' @param freq a vector of allele frequencies
#' @param pool_alleles number of most frequent alleles to pool the frequency for
#' (this is used for calculating modified homozigosity statistics such as H12)
expectedHomozygosity <- function(freq, pool_alleles = 1){

  # check input
  if(!is.numeric(freq)) stop("freq must be numeric vector.")
  if(pool_alleles < 1 | pool_alleles >= length(freq)) stop("pool_alleles not valid.")

  # sort frequencies
  freq <- sort(freq, decreasing = TRUE)

  # get alleles to pool together
  freq1 <- sum(freq[1:pool_alleles])
  freq2 <- freq[-(1:pool_alleles)]

  # expected homozygosity
  hom <- freq1^2 + sum(freq2^2)

  return(hom)
}

#' Calculate several diversity metrics
#' @param x a tabble of frequencies as output by readHarp
estimateDiversity <- function(x){
  x %>%
    group_by(sample, nitrate, selection, rep, chrom, pos, start, end) %>%
    summarise(hap_sum = sum(freq),
              hap_hom1 = expectedHomozygosity(freq, 1),
              hap_hom12 = expectedHomozygosity(freq, 2),
              hap_hom123 = expectedHomozygosity(freq, 3),
              hap_hom1234 = expectedHomozygosity(freq, 4),
              hap_hom2 = hap_hom1 - sort(freq, decreasing = TRUE)[1]^2,
              hap_shannon = -sum(freq[freq > 0] * log(freq[freq > 0])),
              hap_diversity = exp(hap_shannon),
              hap_nalleles = sum(freq > 1/200),
              freqs_csv = paste(acc, freq, sep = ",", collapse = "\n")) %>%
    ungroup()
}


#
# Read HARP data ----
#
# Read each of the window sizes to a list
harp_freqs <- readHarp(paste0(opt$window, "kb.freqs"), opt$indir)

# Calculate diversity metrics
harp_diversity <- estimateDiversity(harp_freqs)


#
# Tidy data ----
#
# Remove resequenced pools from dataset (they are highly correlated)
reseq_samples <- unique(harp_freqs$sample)
reseq_samples <- reseq_samples[grep("_reseq", reseq_samples)]
reseq_samples <- gsub("_reseq", "", reseq_samples)

freqs <- harp_freqs %>%
    filter(!(sample %in% reseq_samples)) %>%
    mutate(rep = str_replace(rep, "_reseq", ""))

diversity <- harp_diversity %>%
    filter(!(sample %in% reseq_samples)) %>%
    mutate(rep = str_replace(rep, "_reseq", ""))


#
# Write files ----
#
# Write window summaries
write_csv(freqs, paste0(opt$outdir, "/haplotype_freq_", opt$window, "kb.csv"))
write_csv(diversity, paste0(opt$outdir, "/haplotype_diversity_", opt$window, "kb.csv"))

