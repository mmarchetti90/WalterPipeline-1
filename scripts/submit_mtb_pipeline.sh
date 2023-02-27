#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --account=owner-guest
#SBATCH --partition=kingspeak-guest
#SBATCH -o mtb-call_run-out-%j
#SBATCH -e mtb-call_run-err-%j
#SBATCH --mail-type=ALL
#SBATCH --mail-user=katharine.walter@hsc.utah.edu
#SBATCH --job-name=test

module load nextflow/20.10 singularity/3.8.7

export WORKDIR=/uufs/chpc.utah.edu/common/home/walter-group1/tb/WalterPipeline
export TMPDIR=/scratch/general/nfs1/u6045141/tmp
export NXF_SINGULARITY_CACHEDIR=$WORKDIR/images

echo "Job started at $(date)"
cd $WORKDIR
nextflow run main.nf -profile singularity
cd $HOME
echo "Job ended at $(date)"
