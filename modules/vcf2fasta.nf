process ConvertVCF {

  // Convert single sample VCF to fasta

  label 'slurm'

  publishDir "${projectDir}/results/${batch}/${sample_id}/fasta", mode: "copy", pattern: "*_${params.aller}.fa"

  input:
  each path(reference)
  tuple val(sample_id), path(unfilt_vcf), val(batch), val(run)

  output:
  path "${sample_id}_${params.caller}.fa"

  """
  # Get sample name for correct genotype
  sample=\$(bcftools query -l ${unfilt_vcf} )

  # Index vcf
  tabix -p vcf ${unfilt_vcf}

  # Consensus with no masking (exclude indels).
  bcftools consensus --include 'TYPE!="indel"' --fasta-ref ${reference} --sample \${sample} --absent 'N' --missing 'N' ${unfilt_vcf} | \
  sed "s/>NC_000962.3 Mycobacterium tuberculosis H37Rv, complete genome/>$sample_name/g" > ${sample_id}_${params.caller}.fa
  """

}