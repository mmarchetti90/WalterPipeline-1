#!/bin/bash

########################################################################################
#### Download references & index reference genomes for Mtb variant calling pipeline ####
########################################################################################

# Define local container tool: docker (default), podman, or singularity.
# N.B. If using Podman, choose "docker"
container=${1:-"docker"}

# Define locations
res_dir=resources/
ref_dir=refs/
bwa_index_dir=bwa_index/
bowtie2_index_dir=bowtie2_index/
gatk_dictionary_dir=gatk_dictionary/
kraken2_db_dir=kraken_db/

# Define main variables
ref_name="H37Rv.fasta"
gatk_dictionary_name="H37Rv.dict"
bowtie_index_prefix="H37Rv"
snpeff_url=https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
ncbi_id="NC_000962.3"

# Make necessary directories
mkdir -p ${res_dir}
cd ${res_dir}
mkdir -p ${ref_dir}
mkdir -p ${bwa_index_dir}
mkdir -p ${bowtie2_index_dir}
mkdir -p ${gatk_dictionary_dir}
mkdir -p ${kraken2_db_dir}

# Define run command and options
if [ "$container" = "docker" ]
then

	run_command="run"
	bind_option="-v $(pwd):/home"
	other_options="--rm"
	image="ksw9/mtb-call"

elif [ "$container" = "podman" ]
then

	run_command="run"
	bind_option="-v $(pwd):/home"
	other_options="--rm"

else

	run_command="exec"
	bind_option=""
	other_options="--cleanenv"
	image="docker://ksw9/mtb-call"
	
fi

# Requires: entrez-direct (conda), bwa, GATK, samtools, kraken2, all present in the container.
# Run each step within the image. Need to mount local directory so that the resources are downloaded locally, not just in the container.

## 1. Reference genome ##
# Download H37Rv reference fasta (alternatively, use efetch)
${container} ${run_command} ${bind_option} ${other_options} ${image} /bin/bash -c "esearch -db nucleotide -query ${ncbi_id} | efetch -format fasta > ${ref_dir}${ref_name}"

# bwa index reference
${container} ${run_command} ${bind_option} ${other_options} ${image} bwa index ${ref_dir}${ref_name}
mv ${ref_dir}*.{amb,ann,bwt,pac,sa} ${bwa_index_dir}

# bowtie2 index reference
${container} ${run_command} ${bind_option} ${other_options} ${image} bowtie2-build ${ref_dir}${ref_name} ${bowtie2_index_dir}${bowtie_index_prefix}

# Create fasta index file
${container} ${run_command} ${bind_option} ${other_options} ${image} samtools faidx ${ref_dir}${ref_name}

# create GATK reference dictionary
${container} ${run_command} ${bind_option} ${other_options} ${image} gatk CreateSequenceDictionary -R ${ref_dir}${ref_name}
mv ${ref_dir}${gatk_dictionary_name} ${gatk_dictionary_dir}

## 2. Masking bed file ##
# PPE masking bed file (H37Rv_ppe.bed.gz) and VCF PPE annotation (ppe_hdr.txt) are in resources/bed hosted on Github

## 3. Annotation information ##
# Download SnpEff for gene annotation.
${container} ${run_command} ${bind_option} ${other_options} ${image} wget ${snpeff_url}

# Unzip file
unzip snpEff_latest_core.zip
rm snpEff_latest_core.zip

# Download the updated M. tuberculosis annotations.
${container} ${run_command} ${bind_option} ${other_options} ${image} java -jar snpEff/snpEff.jar download Mycobacterium_tuberculosis_h37rv

## 4. Kraken2 Database ##
# download kraken2 database (requires ~100G) 
# https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown
${container} ${run_command} ${bind_option} ${other_options} ${image} kraken2-build --standard --threads 24 --db ${kraken2_db_dir}

## 5. Create reads list input with full paths to test data.
cd ../test_reads
touch ../resources/reads_list.tsv
echo "sample\tfastq_1\tfastq_2\tbatch\trun" >> ../resources/reads_list.tsv
echo "PipelineTesting\t$(pwd)/test_R1_001.fastq.gz\t$(pwd)/test_R1_001.fastq.gz\tTest\tTestRun" >> ../resources/reads_list.tsv
