process VariantsLoFreq {

  // Variant calling with LoFreq
  
  label 'slurm'

  publishDir "${projectDir}/results/${batch}/${sample_id}/vars", mode: "copy", pattern: "*{_lofreq_unfilt,_lofreq_filt}.vcf.gz"

  input:
  each path(reference)
  each path(reference_index)
  tuple val(sample_id), path(bam), val(batch)

  output:
  tuple val(sample_id), path("${sample_id}_lofreq_unfilt.vcf.gz"), val(batch), emit: lofreq_unfilt_vcf
  tuple val(sample_id), path("${sample_id}_lofreq_filt.vcf.gz"), val(batch), emit: lofreq_filt_vcf

  """
  # Indexing bam
  samtools index ${bam}

  # Call variants with LoFreq, no filter
  lofreq call-parallel --call-indels --pp-threads \$SLURM_CPUS_ON_NODE --no-default-filter -f ${reference} -o ${sample_id}_lofreq_unfilt.vcf ${bam}

  # Bgzipping
  bgzip ${sample_id}_lofreq_unfilt.vcf

  # Call variants with LoFreq, default filter
  lofreq call-parallel --call-indels --pp-threads \$SLURM_CPUS_ON_NODE -f ${reference} -o ${sample_id}_lofreq_filt.vcf ${bam}

  # Bgzipping
  bgzip ${sample_id}_lofreq_filt.vcf
  """

}