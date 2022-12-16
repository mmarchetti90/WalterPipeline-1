process MapReads_BWA {

  // Map reads and remove duplicates in same step
  
  label 'slurm'

  publishDir "${projectDir}/results/${batch}/${sample_id}/bams", mode: "copy", pattern: "*_merged_rmdup.bam"

  input:
  each path(reference_fasta)
  path(bwa_index)
  tuple val(sample_id), path(read1), path(read2), val(batch), val(run)

  output:
  tuple val(sample_id), path("${sample_id}_merged_rmdup.bam"), val(batch), val(run), emit: bam_files

  """
  # Get machine id_lane from SRA read identifier (old Illumina fastq format)
  seqid=\$(zcat ${read1} | head -n 1)
  seqid="\$(echo \$seqid | cut -d' ' -f1)"
  seqid="\$(echo \$seqid | cut -d':' -f3)"
  id_lane=\${seqid:-readgroup1} # default is "readgroup1"

  # Mapping with BWA
  bwa mem -t 7 ${reference_fasta} ${read1} ${read2} > temp.sam

  # Convert sam to bam 
  sambamba view -t 7 -S -h temp.sam -f bam -o temp.bam

  # Add/replace read groups for post-processing with GATK
  picard AddOrReplaceReadGroups \
  -INPUT temp.bam \
  -OUTPUT temp_rg.bam \
  -RGID \${id_lane} \
  -RGLB library1 \
  -RGPU \${id_lane} \
  -RGPL "illumina" \
  -RGSM ${sample_id}

  # Sort the BAM 
  sambamba sort --tmpdir . temp_rg.bam

  # Index BAM
  sambamba index temp_rg.sorted.bam

  # Remove duplicates. (-r = remove)
  sambamba markdup -r -p -t \$SLURM_CPUS_ON_NODE --tmpdir . temp_rg.sorted.bam ${sample_id}_merged_rmdup.bam
  """

}