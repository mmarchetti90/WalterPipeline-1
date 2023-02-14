#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --account=owner-guest
#SBATCH --partition=redwood-guest
#SBATCH --nodes=1
#SBATCH --mem=0    
#SBATCH --mail-type=ALL
#SBATCH --mail-user=u6045141@utah.edu
#SBATCH --job-name=TestRWWRC

export WORKDIR=/scratch/general/pe-nfs1/u6045141/WalterPipeline
export TMPDIR=$WORKDIR/tmp
export NXF_SINGULARITY_CACHEDIR=$TMPDIR
export SINGULARITY_CACHEDIR=$TMPDIR
export SINGULARITY_TMPDIR=$TMPDIR

module load singularity
module load nextflow

echo "Job started at $(date)"
cd $WORKDIR
nextflow run main.nf -profile singularity
cd $HOME
echo "Job ended at $(date)"
