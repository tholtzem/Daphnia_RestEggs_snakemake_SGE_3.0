# Daphnia_RestEggs_snakemake_SGE_3.0

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.6.0-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)

======================================================

Snakemake workflow for the analysis of DNA from *Daphnia*-resting eggs on the leo4 HPC cluster with the SGE batch job submission system. 

The snakemake workflow was modified from the [Daphnia_RestEggs_snakemake_pbs_2.0](https://github.com/tholtzem/Daphnia_RestEggs_snakemake_pbs_2.0), which was initially based on the [ta_dna_snakemake_pbs Tutorial](https://github.com/schimar/ta_dna_snakemake_pbs),the [Snakemake Cluster Tutorial](https://github.com/SchlossLab/snakemake_cluster_tutorial.git) and the [Software Carpentry lesson repository](https://hpc-carpentry.github.io/hpc-python/17-cluster/). For more information on snakemake itself (https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html).

The analysis is based on the [physalia-lcwgs: the Physalia course on Population genomic inference from low-coverage whole-genome sequencing data, Oct 19-22 2020](https://github.com/nt246/physalia-lcwgs).

This workflow starts with mapping to the new reference assembly of [*D. galeata*](https://doi.org/10.1093/gbe/evab267). For details on quality filtering and cleaning reads from non-target DNA with Kraken2 see [Daphnia_RestEggs_snakemake_pbs_2.0](https://github.com/tholtzem/Daphnia_RestEggs_snakemake_pbs_2.0). Please note that I used PBS submission scripts in the previous workflow.

======================================================

## conda/mamba and other [dependencies](https://github.com/schimar/ta_dna_snakemake_pbs/blob/main/envs/s21.yaml)   

Mamba (https://github.com/mamba-org/mamba) is a reimplementation of the conda package manager in C++.

```
# Load Anaconda on cluster (here mach2):
module load Anaconda3/2021.04/miniconda-base-2021.04

# To use conda commands in your current shell session, first do:
source $UIBK_CONDA_PROFILE

# Create environment from yaml file (in envs/):
conda init bash
mamba env create -f envs/s21.yaml

# Activate the environment
conda activate eggs3.0

# if you've added new software to install to the conda environment, then you can update:
mamba env update --name eggs3.0 --file envs/s21.yaml

# In case you want to remove the conda environment
conda env remove -n eggs3.0

```

