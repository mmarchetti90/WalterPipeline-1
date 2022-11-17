process MapReads_Bowtie {

  // Map reads and remove duplicates in same step
  
  label 'slurm'

  publishDir "${projectDir}/results/${batch}/${sample_id}/bams", mode: "copy", pattern: "*_merged_rmdup.bam"

  input:
  each path bowtie_index
  tuple val(sample_id), path(read1), path(reads2), val(batch), val(run)

  output:
  tuple val(sample_id), path("${sample_id}_merged_rmdup.bam"), val(batch), val(run) emit: bam_files

  """
  # Get machine id_lane from SRA read identifier (old Illumina fastq format)
  seqid=\$(zcat ${read1} | head -n 1)
  seqid="\$(echo $seqid | cut -d' ' -f1)"
  seqid="\$(echo $seqid | cut -d':' -f3)"
  id_lane=\${seqid:-readgroup1} # default is "readgroup1"

  # Alignment with Bowtie2
  bowtie2 --threads \$SLURM_CPUS_ON_NODE -X 1100 -x ${params.bowtie_index_prefix} -1 ${read1} -2 ${read2} -S temp.sam

  # Convert sam to bam 
  sambamba view -t \$SLURM_CPUS_ON_NODE -S -h temp.sam -f bam -o temp.bam

  # Add/replace read groups for post-processing with GATK
  picard AddOrReplaceReadGroups \
  INPUT=temp.bam \
  OUTPUT=temp_rg.bam \
  RGID=\${id_lane} \
  RGLB=library1 \
  RGPU=\${id_lane} \
  RGPL="illumina" \
  RGSM=${sample_id}

  # Sort the BAM 
  sambamba sort temp_rg.bam

  # Index BAM
  sambamba index temp_sorted.bam

  # Remove duplicates. (-r = remove)
  sambamba markdup -r -p -t \$SLURM_CPUS_ON_NODE temp_sorted.bam ${sample_id}_merged_rmdup.bam
  """

}