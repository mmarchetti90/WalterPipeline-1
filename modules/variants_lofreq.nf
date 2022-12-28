process VariantsLoFreq {

  // Variant calling with LoFreq
  
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
  lofreq call-parallel --call-indels --pp-threads \$SLURM_CPUS_ON_NODE --no-default-filter -f ${reference} -o ${sample_id}_lofreq.vcf ${bam}

  # Gzipping
  bgzip ${sample_id}_lofreq.vcf
  """

}
