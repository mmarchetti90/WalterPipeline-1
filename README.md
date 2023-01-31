# *M. tuberculosis*  variant identification pipeline

Pipeline for *M. tuberculosis* variant identification from short-read data for epidemiology and phylogenetics. Briefly, this pipeline takes raw short-read Illumina data, pre-processes reads (read adapter trimming and taxonomic filtering), maps reads to a provided reference genome, calls variants, and outputs consensus FASTA sequences for downstream applications. The pipeline is written in Nextflow, so modules are adaptable. User options allow for tailoring of the pipeline such as setting custom filters and choosing a reference genome.

## Installation & set-up

1. Clone Github repo.
```
git clone https://github.com/ksw9/WalterPipeline.git
```

2. Load your HPC's container tool (i.e. Docker or Singularity) and nextflow. (Some clusters may have these pre-loaded.)
```
module load singularity # or docker
module load nextflow
```

3. Run script to download references and resources, specify docker/sigularity (these, especially the Kraken2 database, are too large to include elsewhere). The Kraken2 database requires ~100G of space; users with more limited memory might consider a different database.
```
# Run download_refs.sh.
cd WalterPipeline 
./scripts/download_refs.sh singularity # or docker
```
This should populate your resources directory with all required references and databases.

3. Modify the config file (nextflow.config).
  - update resources_dir (full path to directory resources)
  - update clusterOptions parameter to make specific to cluster.

## Usage
1. Run the pipeline on the test data (truncated FASTQ files) included in the test_data directory. Include any user options here. The Docker image will be pulled automatically by running the pipeline the first time.
```
nextflow run main.nf
```

2. Run the pipeline on user data. 
  - Create a tab-delimited file with sample name, full path to FASTQ read 1, full path to FASTQ read 2, batch name, run name (format like data/reads_list.tsv). 
  - Update the nextflow.config so that the reads_list parameter is now defined by the new list. 
  - Run the pipeline.
```
nextflow run main.nf
```

## Options

There are several user options which can be modified on the command line or in the nextflow.config file (command line options take precedence).
- mapper (bwa/bowtie2): defines mapping algorithm to be used (default = bwa).
- run_lofreq (true/false): In addition to calling variants with GATK, will call low frequency minority variants with LoFreq.
- depth_threshold: defines minimum site depth for calling an allele (either variant or reference) that will be applied to generate a consensus sequence (default = 5)
- qual_threshold: defines the minimum site quality score for calling an allele (either variant or reference) that will be applied to generate a consensus sequence (default = 20)
- ploidy: defines ploidy for GATK variant calling (currently, only tested for ploidy = 1)
- threads: defines available threads (default = 4)
- nextseq (true/false): Use of NextSeq sequencing platform? (default = false). Nextseq has been found to [overcall](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md) G bases at the 3' end; if this option is turned on, TrimGalore will ignore quality scores of G bases in the trimming step. 
- nextseq_qual_threshold: If the above parameter is true, defines the quality threshold for trimming (default = 20).

## Outputs

All outputs are stored in the results directory, within the project directory. Directory structure mirrors the input reads file, with directories organized by sequencing run, then sample.
```
├── results
│   ├── test_data/test (example organized by sequencing batch, then sample) 
|   │   ├──trim
|   │   ├──kraken
|   │   ├──bams
|   │   ├──vars
|   │   ├──fasta
|   │   ├──stats
```
## Example data

- Truncated paired-end fastq files are in the test_data directory.
- An input sample .tsv file list is located at config/test_data.tsv.

## Troubleshooting

- Singularity uses the $HOME directory as the default cache. Set to TMPDIR via: 
``` 
export SINGULARITY_CACHEDIR=$TMPDIR
export SINGULARITY_TMPDIR=$TMPDIR
```
