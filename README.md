# *M. tuberculosis*  variant identification pipeline

Pipeline for *M. tuberculosis* variant identification from short-read data for epidemiology and phylogenetics. Briefly, this pipeline takes raw short-read Illumina data, pre-processes reads (read adapter trimming and taxonomic filtering), maps reads to a provided reference genome, calls variants, and outputs consensus FASTA sequences for downstream applications. The pipeline is written in Nextflow, so modules are adaptable. User options allow for tailoring of the pipeline such as setting custom filters and choosing a reference genome.

## Installation & set-up

1. Download Docker image. This is a containerized pipeline, so the Docker image contains the software and versions required for analysis.
```
docker pull YYY 
```

2. Run script to download references and resources (these, especially the Kraken2 database, are too large to include elsewhere). The Kraken2 database requires ~100G of space; users with more limited memory might consider a different database.
```
# Clone Github (includes scripts and small, pipeline-specific resources).
git clone https://github.com/ksw9/WalterPipeline.git

# Run download_refs.sh.
cd WalterPipeline (update w/name of pipeline)
./scripts/download_refs.sh
```
This should populate your resources directory with all required references and databases.

3. Modify the config file (nextflow.config).
  - update the path to the Docker image
  - update all paths to resources

## Usage
1. Run the pipeline on the test data (truncated FASTQ files) included in the test_data directory. Include any user options here. 
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



