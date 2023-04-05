#!/bin/bash
#SBATCH -A LEYSER-SL2-CPU
#SBATCH -J simulations
#SBATCH -D /rds/project/ol235/rds-ol235-leyser-hpc/projects/2020_Tavares_NitrateSelection/supplementary_data/
#SBATCH -o logs/simulations.log
#SBATCH -p skylake
#SBATCH -c 1
#SBATCH --mem-per-cpu=1000MB
#SBATCH -t 00:50:00
#SBATCH -a 1-100  # replicate simulations

# This should ideally be incorporated into snakemake...
source $(conda info --base)/etc/profile.d/conda.sh
conda activate simupop

# loop through parameters
# while read -r PARAMS
# do
#   python workflow/scripts/simulations/single_replicate_simulation.py \
#     --n_selected_loci $(echo $PARAMS | cut -d "," -f 1) \
#     --selected_effect $(echo $PARAMS | cut -d "," -f 2) \
#     --n_adv_alleles $(echo $PARAMS | cut -d "," -f 3) \
#     --outdir "data/intermediate/simulations" \
#     --seed "20200916$SLURM_ARRAY_TASK_ID"
# done < "data/intermediate/simulations/parameter_list.csv"

for NSEL in 1 3 5 10 20 30 40 50 60
do
  for EFF in 0.05 0.1 0.2 0.3 0.5 0.7 1
  do
    for NADV in 1 2 4 8
    do
      # echo "$NSEL-$EFF-$NADV"
      python workflow/scripts/simulations/single_replicate_simulation.py \
        --n_selected_loci $NSEL \
        --selected_effect $EFF \
        --n_adv_alleles $NADV \
        --outdir "data/intermediate/simulations/" \
        --suffix "${NSEL}-${EFF}-${NADV}-seed${SLURM_ARRAY_TASK_ID}" \
        --seed "20200916$SLURM_ARRAY_TASK_ID"
        
    done
  done
done

# approximately infinitesimal - 500 loci, i.e. 100 per chromosome
NSEL="500"
for EFF in 0.01 0.03 0.05 0.08 0.1
do
  for NADV in 1 2 4 8
  do
    # echo "$NSEL-$EFF-$NADV"
    python workflow/scripts/simulations/single_replicate_simulation.py \
      --n_selected_loci 500 \
      --selected_effect $EFF \
      --n_adv_alleles $NADV \
      --outdir "data/intermediate/simulations/" \
      --suffix "${NSEL}-${EFF}-${NADV}-seed${SLURM_ARRAY_TASK_ID}" \
      --seed "20200916$SLURM_ARRAY_TASK_ID"
  done
done

# # 5 loci clumped into one chromosome
# for EFF in 0.1 0.2 0.3 0.5
# do
#   for NADV in 1 2 4 8
#   do
#     # echo "$NSEL-$EFF-$NADV"
#     python workflow/scripts/simulations/single_replicate_simulation.py \
#       --loc_selected_loci "51,101,152,202,253" \
#       --selected_effect $EFF \
#       --n_adv_alleles $NADV \
#       --outdir "data/intermediate/simulations/" \
#       --suffix "5in1-${EFF}-${NADV}-seed${SLURM_ARRAY_TASK_ID}" \
#       --seed "20200916$SLURM_ARRAY_TASK_ID"
      
#   done
# done

# # 10 loci clumped into two chromosomes
# for EFF in 0.1 0.2 0.3 0.5
# do
#   for NADV in 1 2 4 8
#   do
#     # echo "$NSEL-$EFF-$NADV"
#     python workflow/scripts/simulations/single_replicate_simulation.py \
#       --loc_selected_loci "51,101,152,202,253,354,404,455,505,556" \
#       --selected_effect $EFF \
#       --n_adv_alleles $NADV \
#       --outdir "data/intermediate/simulations/" \
#       --suffix "10in2-${EFF}-${NADV}-seed${SLURM_ARRAY_TASK_ID}" \
#       --seed "20200916$SLURM_ARRAY_TASK_ID"
#   done
# done

# # 5 loci closely linked in one chromosome
# for EFF in 0.1 0.2 0.3 0.5
# do
#   for NADV in 1 2 4 8
#   do
#     # echo "$NSEL-$EFF-$NADV"
#     python workflow/scripts/simulations/single_replicate_simulation.py \
#       --loc_selected_loci "26,51,76,101,127" \
#       --selected_effect $EFF \
#       --n_adv_alleles $NADV \
#       --outdir "data/intermediate/simulations/" \
#       --suffix "5in1linked-${EFF}-${NADV}-seed${SLURM_ARRAY_TASK_ID}" \
#       --seed "20200916$SLURM_ARRAY_TASK_ID"
#   done
# done

# # 5 loci in one chromosome and 1 in each of the others
# for EFF in 0.1 0.2 0.3 0.5
# do
#   for NADV in 1 2 4 8
#   do
#     # echo "$NSEL-$EFF-$NADV"
#     python workflow/scripts/simulations/single_replicate_simulation.py \
#       --loc_selected_loci "51,101,152,202,253,455,650,883,1067" \
#       --selected_effect $EFF \
#       --n_adv_alleles $NADV \
#       --outdir "data/intermediate/simulations/" \
#       --suffix "5in1plus4-${EFF}-${NADV}-seed${SLURM_ARRAY_TASK_ID}" \
#       --seed "20200916$SLURM_ARRAY_TASK_ID"
#   done
# done
