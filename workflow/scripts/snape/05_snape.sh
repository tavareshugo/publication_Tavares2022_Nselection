#!/bin/bash
#SBATCH --job-name=snape
#SBATCH --workdir=/home/hugot/projects/20150501_selection_experiment/poolseq/
#SBATCH -N 1
#SBATCH --ntasks 48
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G
#SBATCH -o ./scripts/shell/05_snape.stdout


# sbatch /home/hugot/projects/20150501_selection_experiment/poolseq/scripts/shell/05_snape.sh


###
# Environment variables ####
###

# List of pileup files
PILEUPS=$(find ./mapping/merged -type f -name "*.pileup")

# Make output directory if non-existent
mkdir -p ./diversity_analysis/allele_frequency/all_snps



###
# Run snape-pooled for each PILEUP
###

# Note: the priors for SNAPE were determined empirically:
## median theta ~ 0.004 across all populations, calculated from popoolation software
## divergence calculated from rounded number of SNPs per bp using the list of known SNPs

# Temporary script to run SNAPE in parallel with srun
echo '#!/bin/bash
PILEUP=$1
SAMPLE=$(basename $PILEUP .pileup | sed "s/\.MQ20\.BQ30\.DP10to300//g")

cat $PILEUP | \
snape-pooled -nchr 400 -theta 0.004 -D 0.02 -fold unfolded -priortype informative | \
grep -v "\*" > ./diversity_analysis/allele_frequency/all_snps/${SAMPLE}.snape.tsv
' > ./diversity_analysis/allele_frequency/run_snape.sh
chmod 775 ./diversity_analysis/allele_frequency/run_snape.sh


for PILEUP in $PILEUPS
do
  
  # Extract sample name from pileup file name
  SAMPLE=$(basename $PILEUP .pileup | sed 's/\.MQ20\.BQ30\.DP10to300//g')
  
  # Launch snape in parallel using srun
  srun -n 1 --exclusive -o ./diversity_analysis/allele_frequency/${SAMPLE}.log \
    ./diversity_analysis/allele_frequency/run_snape.sh $PILEUP &
  
done
wait

# remove temporary script
rm ./diversity_analysis/allele_frequency/run_snape.sh
