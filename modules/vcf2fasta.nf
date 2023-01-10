process ConvertVCF {

  // Convert single sample VCF to fasta

  label 'slurm'

  publishDir "${projectDir}/results/${batch}/${sample_id}/fasta", mode: "copy", pattern: "*_${params.variant_caller}.fa"

  input:
  each path(reference)
  tuple val(sample_id), path(vcf), val(batch), val(run)

  output:
  path "${sample_id}_${params.variant_caller}.fa"

  """
  # Index vcf
  tabix -p vcf ${vcf}

  # Consensus with no masking (exclude indels).
  # N.B. The vcf files come from individual samples, so no need to specify --sample in bcftools consensus (also, LoFreq does not store sample name info in the vcf).
  bcftools consensus --include 'TYPE!="indel"' --fasta-ref ${reference} --missing 'N' --absent 'N' ${vcf} | sed "s/>NC_000962.3 Mycobacterium tuberculosis H37Rv, complete genome/>${sample_id}/g" > ${sample_id}_${params.variant_caller}.fa
  """

}