#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --account=[add account]
#SBATCH --partition=[add partition]
#SBATCH --nodes=1
#SBATCH --mem=0    
#SBATCH --mail-type=ALL
#SBATCH --mail-user=[add email]
#SBATCH --job-name=test

export WORKDIR=[add working directory]
export TMPDIR=[add temp directory]
export NXF_SINGULARITY_CACHEDIR=$WORKDIR/images
export SINGULARITY_CACHEDIR=$WORKDIR/images
export SINGULARITY_TMPDIR=$TMPDIR

module load singularity
module load java nextflow

echo "Job started at $(date)"
cd $WORKDIR
nextflow run main.nf -profile singularity
cd $HOME
echo "Job ended at $(date)"
