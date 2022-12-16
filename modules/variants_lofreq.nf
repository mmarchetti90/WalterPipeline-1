process VariantsLoFreq {

  // Run AMR prediction tool. Removing for now
  
  label 'slurm'

  publishDir "${projectDir}/results/${batch}/${sample_id}/vars", mode: "copy", pattern: "*_lofreq.vcf.gz"

  input:
  each path(reference)
  each path(reference_index)
  tuple val(sample_id), path(bam), val(batch), val(run)

  output:
  tuple val(sample_id), path("${sample_id}_lofreq.vcf.gz"), val(batch), val(run), emit: lofreq_vcf

  """
  # Indexing bam
  samtools index ${bam}

  # Call variants with LoFreq
  lofreq call-parallel --pp-threads \$SLURM_CPUS_ON_NODE -f ${reference} -o ${sample_id}_lofreq.vcf ${bam}

  # Gzipping
  bgzip ${sample_id}_lofreq.vcf
  """

}