# *M. tuberculosis*  variant identification pipeline

Pipeline for *M. tuberculosis* variant identification from short-read data for epidemiology and phylogenetics. Briefly, this pipeline takes raw short-read Illumina data, pre-processes reads (read adapter trimming and taxonomic filtering), maps reads to a provided reference genome, calls variants, and outputs consensus FASTA sequences for downstream applications. The pipeline is written in Nextflow, so modules are adaptable. User options allow for tailoring of the pipeline such as setting custom filters and choosing a reference genome.

## Installation & set-up

1. Download Docker image. This is a containerized pipeline, so the Docker image contains the software and versions required for analysis.
```
docker pull YYY 
```

2. Run script to download references and resources (these, especially the Kraken2 database, are too large to include elsewhere). The Kraken2 database requires ~100G of space; users with more limited memory might consider a different database.
```
# Download script

#
chmod 755 script

./script
```

3. Modify the config file (nextflow.config) and sample input sheet (xxx). 


## Usage
1. Run the pipeline. Include any user options here. 
```
nextflow run main.nf
```
