#
# Fig S14
#

library(seqinr)
library(magrittr)
library(ape)
library(pegas)


max2 <- read.alignment("./data/raw/max2_mutation/MAX2_sequence.fas",
                       format = "fasta") %>%
  as.DNAbin()

max2_hap <- max2 %>% haplotype()
max2_hap <- sort(max2_hap, what = "labels")
max2_hapnet <- max2_hap %>% haploNet()

# Function to check which accessions have each haplotype
# see: https://stackoverflow.com/questions/25755930/how-to-plot-pie-charts-in-haplonet-haplotype-networks-pegas/25756818#25756818
countHap <- function(hap = h, dna = x){
  with(
    stack(setNames(attr(hap, "index"), rownames(hap))),
    table(hap = ind, pop = attr(dna, "dimnames")[[1]][values])
  )
}

labs <- countHap(max2_hap, max2) %>%
  apply(1, function(i){
    paste(names(i)[which(i == 1)], collapse = "\n")
  })

attr(max2_hapnet, "labels") <- labs

pdf("./figures/S14Fig.pdf", width = 7.5, height = 4)
plot(max2_hapnet,
     size = attr(max2_hapnet, "freq") * 1.5,
     fast = TRUE, labels = TRUE,
     font = 1, cex = 0.9, threshold = 0,
     show.mutation = 3,
     bg = "lightgrey")
dev.off()
