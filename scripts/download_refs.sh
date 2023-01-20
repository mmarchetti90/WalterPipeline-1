#!/bin/bash

########################################################################################
#### Download references & index reference genomes for Mtb variant calling pipeline ####
########################################################################################

# Define local container tool (default is docker).
container=${1:-"docker"}   

# Define locations
ref_dir=resources/
mkdir -p ${ref_dir}
cd ${ref_dir}
ref_path=${ref_dir}/refs/H37Rv.fasta
ncbi_id="NC_000962.3"
image="ksw9/mtb-call"

# Requires: entrez-direct (conda), bwa, GATK, samtools, kraken2, all present in the container.
# Run each step within the image. Need to mount local directory so that the resources are downloaded locally, not just in the container.
## 1. Reference genome ##
# Download H37Rv reference fasta (alternatively, use efetch)
${container} run -v $(pwd)/${ref_dir} esearch -db nucleotide -query ${ncbi_id} | efetch -format fasta > ${ref_path}

# bwa index reference
${container} run -v $(pwd)/${ref_dir} bwa index ${ref_path}

# Create fasta index file
${container} run -v $(pwd)/${ref_dir}  samtools faidx ${ref_path}

# create GATK reference dictionary
${container} run -v $(pwd)/${ref_dir} gatk CreateSequenceDictionary -R ${ref_path}

## 2. Masking bed file ##
# PPE masking bed file (H37Rv_ppe.bed.gz) and VCF PPE annotation (ppe_hdr.txt) are in resources/bed hosted on Github

## 3. Annotation information ##
# Download SnpEff for gene annotation.
cd ${ref_dir}/snpEff
${container} run -v $(pwd)/${ref_dir} wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip

# Unzip file
unzip snpEff_latest_core.zip

# Download the updated M. tuberculosis annotations.
${container} run -v $(pwd)/${ref_dir} java -jar snpEff.jar download Mycobacterium_tuberculosis_h37rv

## 4. Kraken2 Database ##
# download kraken2 database (requires ~100G) 
# https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown
${container} run -v $(pwd)/${ref_dir} kraken2-build --standard --threads 24 --db $DBNAME

## 5. Create reads list input with full paths to test data. (To update!)

