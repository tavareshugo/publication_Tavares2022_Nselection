library(seqinr)

fas <- read.fasta("data/external/branch_regulators.fas")

# keep only major isoform
fas <- fas[str_detect(names(fas), ".1$")]

# function to count possible mutations resulting in stop codon
count_stop_mutation <- function(codon){
  codon <- toupper(codon)

  if(codon %in% c("TAG", "TAA", "TGA")) return(0)

  n <- 0

  # check for possible mutations that change a codon to a stop
  if(str_detect(codon, "^TA")) n <- n + 2
  if(str_detect(codon, "^TG")) n <- n + 1
  if(str_detect(codon, "AG$")) n <- n + 1
  if(str_detect(codon, "GA$")) n <- n + 1
  if(str_detect(codon, "AA$")) n <- n + 1
  if(str_detect(codon, "T[A,C,G,T]G")) n <- n + 1
  if(str_detect(codon, "T[A,C,G,T]A")) n <- n + 2

  return(n)
}

lapply(fas, function(gene){
  nonsyn_codons <- map_dbl(splitseq(gene), count_stop_mutation)
}) %>%
  unlist() %>%
  sum()
