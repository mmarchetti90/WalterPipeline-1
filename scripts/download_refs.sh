#!/bin/bash

########################################################################################
#### Download references & index reference genomes for Mtb variant calling pipeline ####
########################################################################################

# Requires: entrez-direct (conda), bwa, GATK, samtools, kraken2

# Define locations
ref_dir=resources/
mkdir -p ${ref_dir}
cd ${ref_dir}
ref_path=${ref_dir}/refs/H37Rv.fasta
ncbi_id="NC_000962.3"

## 1. Reference genome ##
# Download H37Rv reference fasta (alternatively, use efetch)
esearch -db nucleotide -query ${ncbi_id} | efetch -format fasta > ${ref_path}

# bwa index reference
bwa index ${ref_path}

# Create fasta index file
samtools faidx ${ref_path}

# create GATK reference dictionary
gatk CreateSequenceDictionary -R ${ref_path}

## 2. Masking bed file ##
# PPE masking bed file (H37Rv_ppe.bed.gz) and VCF PPE annotation (ppe_hdr.txt) are in resources/bed hosted on Github

## 3. Annotation information ##
# Download SnpEff for gene annotation.
cd ${ref_dir}/snpEff
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip

# Unzip file
unzip snpEff_latest_core.zip

# Download the updated M. tuberculosis annotations.
java -jar snpEff.jar download Mycobacterium_tuberculosis_h37rv

## 4. Kraken2 Database ##
# download kraken2 database (requires ~100G) 
# https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown
kraken2-build --standard --threads 24 --db $DBNAME