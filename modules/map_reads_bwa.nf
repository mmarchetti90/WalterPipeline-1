process MapReads_BWA {

  // Map reads and remove duplicates in same step
  
  label 'slurm'

  publishDir "${projectDir}/results/${batch}/${sample_id}/bams", mode: "copy", pattern: "*_merged_rmdup.bam"
  publishDir "${projectDir}/results/${batch}/${sample_id}/stats", mode: "copy", pattern: "**_marked_dup_metrics.txt"  

  input:
  each path(reference_fasta)
  path(bwa_index)
  tuple val(sample_id), path(read1), path(read2), val(batch), val(run)

  output:
  tuple val(sample_id), path("${sample_id}_merged_mrkdup.bam"), val(batch), val(run), emit: bam_files
  tuple val(sample_id), path("${sample_id}_marked_dup_metrics.txt"), val(batch), val(run), emit: dup_metrics
  
  """
  # Get machine id_lane from SRA read identifier (old Illumina fastq format)
  seqid=\$(zcat ${read1} | head -n 1)
  seqid="\$(echo \$seqid | cut -d' ' -f1)"
  seqid="\$(echo \$seqid | cut -d':' -f3)"
  id_lane=\${seqid:-readgroup1} # default is "readgroup1"
  
  # Mapping with BWA
  bwa mem -t \$SLURM_CPUS_ON_NODE ${reference_fasta} ${read1} ${read2} > temp.sam

  # Sort and convert to bam
  gatk SortSam \
  I=temp.sam \
  O=temp.bam \
  SORT_ORDER=coordinate
  
  # Convert sam to bam 
  # sambamba view -t \$SLURM_CPUS_ON_NODE -S -h temp.sam -f bam -o temp.bam

  # Add/replace read groups for post-processing with GATK
  picard AddOrReplaceReadGroups \
  -INPUT temp.bam \
  -OUTPUT temp_rg.bam \
  -RGID \${id_lane} \
  -RGLB ${params.library} \
  -RGPU \${id_lane} \
  -RGPL ${params.seq_platform} \
  -RGSM ${sample_id}
  
  # Mark duplicates & produce library complexity metrics. 
  gatk MarkDuplicates \
  I=temp_rg.bam \
  O=${sample_id}_merged_mrkdup.bam  \
  M=${sample_id}_marked_dup_metrics.txt  
      
# Mark duplicates with Spark (also sorts BAM) & produce library complexity metrics. 
# Need to use Java 7 or Java 11 (https://gatk.broadinstitute.org/hc/en-us/community/posts/4417665825307-java-lang-reflect-InaccessibleObjectException-Unable-to-make-field-transient-java-lang-Object-java-util-ArrayList-elementData-accessible-module-java-base-does-not-opens-java-util-to-unnamed-module-3bf44630)
#   gatk MarkDuplicatesSpark \
#   -I temp_rg.bam \
#   -O ${sample_id}_merged_mrkdup.bam \
#   -M ${sample_id}_marked_dup_metrics.txt  
  
  # Base (Quality Score) Recalibration: not done because no 'known variants.'
  
  """

}