#!/bin/bash
#SBATCH -A LEYSER-SL2-CPU
#SBATCH -J simulations
#SBATCH -D /rds/project/ol235/rds-ol235-leyser-hpc/projects/2020_Tavares_NitrateSelection/supplementary_data/
#SBATCH -o logs/simulations.log
#SBATCH -p skylake
#SBATCH -c 1
#SBATCH --mem-per-cpu=1000MB
#SBATCH -t 00:15:00
#SBATCH -a 1-100  # replicate simulations

# This should ideally be incorporated into snakemake...

source activate simupop

# loop through parameters
while read -r PARAMS
do
  python workflow/scripts/simulations/single_replicate_simulation.py \
    --n_selected_loci $(echo $PARAMS | cut -d "," -f 1) \
    --selected_effect $(echo $PARAMS | cut -d "," -f 2) \
    --n_adv_alleles $(echo $PARAMS | cut -d "," -f 3) \
    --outdir "data/intermediate/simulations" \
    --seed "20200916$SLURM_ARRAY_TASK_ID"
done < "data/intermediate/simulations/parameter_list.csv"
