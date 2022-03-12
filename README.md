# Code for Tavares, Readshaw et al

This is still under review and subject to change. 

We're working to ensure consistent file structure for reproducing the analysis and figures in our paper. 

## WGS analysis

Analysis scripts to process Pool-seq data is provided as a _snakemake_ workflow.

To run on a local machine, run the following from the project's directory:

```shell
snakemake --use-conda
```

We ran this pipeline on the University of Cambridge HPC server, as:

```shell
snakemake --use-conda --jobs 5 \
 --cluster "sbatch -A LEYSER-SL2-CPU -J {rulename} -o logs/slurm/{rulename}-job#%j.log \
 -p skylake -c {threads} --mem-per-cpu=5980MB -t {params.runtime}"
```
