#!/bin/bash
#SBATCH --job-name=snape_filter
#SBATCH --workdir=/home/hugot/projects/20150501_selection_experiment/poolseq/
#SBATCH -N 1
#SBATCH --ntasks 24
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH -o ./scripts/shell/06_snape_filtering.stdout


# sbatch /home/hugot/projects/20150501_selection_experiment/poolseq/scripts/shell/06_snape_filtering.sh


###
# Environment variables ####
###

# Working directory
WORKDIR="/home/hugot/projects/20150501_selection_experiment/poolseq/"

# List of snape-pooled files
INPUTS=$(ls ./diversity_analysis/allele_frequency/all_snps/*.snape.tsv)

# Make output directory
mkdir -p ./diversity_analysis/allele_frequency/founder_snps

# create list of known SNPs where first column is chromosome and second column is position
cut -f 1,4 ./starting_population/founder_accession_data/founders_biallelic.map > founder_snps.tsv

for INPUT in $INPUTS
do
  
  # Output file name (same as input but to a different folder)
  OUTPUT=$(echo ${INPUT} | sed 's/all_snps/founder_snps/' | sed 's/\.tsv$/\.csv/')
  
  # Run the snape filtering script to retain only known SNPs
  srun -N 1 -n 1 --exclusive --mem-per-cpu=20G \
    Rscript ~/code/bioPipelines/snape-pooled_parser.R \
      --snape $INPUT \
      --out $OUTPUT \
      --snp_list founder_snps.tsv &
done
wait
