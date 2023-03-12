##################################################
# Prepare data with start population frequencies #
##################################################

# install valr (cannot do it through conda)
if(!("valr" %in% rownames(installed.packages()))){
  install.packages("valr", repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(valr)
  library(readxl)
})


#
# Get MAGIC selection founders ----
#
# get experimental pedigree
ped <- read_csv("./data/raw/phenotypes/compiled_phenotypes_data_clean.csv",
                col_types = "iccccDiiinninccccccclnDiiiDnnncccllc")

# filter to obtain MAGIC parentals
magics <- ped %>%
  # get the first generation individuals
  filter(generation == 0) %>%
  distinct(nitrate_grown, selection, replicate, mother, father) %>%
  pivot_longer(cols = c(mother, father), names_to = "parent", values_to = "id")

# get public names to match with HSRIL names
# Supplementary file 3 https://doi.org/10.1534/genetics.114.170746
# even though the file extension is ".xls" I had to use read_xlsx function
# see https://github.com/tidyverse/readxl/issues/598
magics <- read_xlsx("./data/external/founder_genotypes/genetics.114.170746-3.xls",
           sheet = "seed area", range = "H1:I824") %>%
  # tidy names up
  rename_with(~ str_to_lower(str_replace(., " ", "_"))) %>%
  # join with pedigree table
  right_join(magics, by = c("our_name" = "id"))


#
# Get MAGIC line mosaics ----
#
# TableS2 From Imprialou et al 2017 (https://doi.org/10.1534/genetics.116.192823)
mosaic <- read_table("./data/external/founder_genotypes/imprialou_mosaic.txt")

# Tidy the mosaic table
mosaic <- mosaic %>%
  # rename variables
  rename(chrom = chr, start = from.bp, end = to.bp, accession = acc) %>%
  # remove accession number for simplicity
  mutate(accession = str_replace(accession, "-.*", "")) %>%
  # retain only MAGICs in our dataset
  filter(magic %in% magics$public_name)


#
# Get mosaic overlaps across MAGIC lines ----
#
# Split table by MAGIC line
mosaic_split <- mosaic %>%
  select(magic, chrom, start, end) %>%
  split(.$magic)

# Make recursive intersects to find common intervals
mosaic_intervals <- reduce(mosaic_split, function(i, j){
  intervals <- bed_intersect(i, j) %>%
    # get interval borders
    mutate(start = pmax(start.x, start.y), end = pmin(end.x, end.y)) %>%
    select(chrom, start, end)
})

# Intersect these with the main table
# to get the original mosaic file but only with intervals common to all MAGICs
mosaic_filtered <- bed_intersect(mosaic, mosaic_intervals) %>%
  rename(start = start.y, end = end.y) %>%
  select(-start.x, -end.x) %>%
  rename_with(~ str_replace(., "\\.x", ""))

# add both nomenclature for MAGIC ids
mosaic_filtered <- magics %>%
  distinct(our_name, public_name) %>%
  right_join(mosaic_filtered, by = c("public_name" = "magic")) %>%
  # retain and rename columns of interest
  select(hsril = our_name, magic = public_name,
         chrom, start, end, accession)

# Save results ------------------------------------------------------------

# save tidy mosaic table
mosaic_filtered %>%
  arrange(as.numeric(str_remove(magic, "MAGIC.")), chrom, start) %>%
  write_csv("data/external/founder_genotypes/magic_mosaics.csv")


# Checks ------------------------------------------------------------------

# Confirm that these intervals contain all 324 accessions with mosaics available
mosaic_filtered %>%
  distinct(magic, chrom, accession, start, end) %>%
  bed_cluster() %>%
  count(.id) %>%
  count(n)

# distance between mosaic intervals
# mosaic_filtered %>%
#   distinct(chrom, start, end) %>%
#   group_by(chrom) %>%
#   arrange(start, end) %>%
#   mutate(diff = lead(start) - end) %>%
#   ungroup() %>%
#   select(chrom, start, end, diff) %>%
#   arrange(chrom, start, end) %>%
#   ggplot(aes(diff)) +
#   geom_histogram() +
#   scale_x_log10() +
#   annotation_logticks(sides = "b")

