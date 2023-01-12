process ConvertVCF {

  // Convert single sample VCF to fasta

  label 'slurm'

  publishDir "${projectDir}/results/${batch}/${sample_id}/fasta", mode: "copy"

  input:
  each path(reference)
  tuple val(sample_id), path(vcf), val(batch), val(run)
  each path(bed_path)

  output:
  path "${sample_id}_${params.variant_caller}_mask.fa", emit: masked_fasta
  path "${sample_id}_${params.variant_caller}.fa", emit: unmasked_fasta

  """
  # Index vcf
  tabix -p vcf ${vcf}

  # Consensus with masking of bed_file and by depth/qual thresholds; exclude indels.
  # N.B. The vcf files come from individual samples, so no need to specify --sample in bcftools consensus (also, LoFreq does not store sample name info in the vcf).
  bcftools consensus --mask ${bed_path} --include "(TYPE!='indel' & INFO/DP >= ${params.depth_threshold}) & (QUAL >= ${params.qual_threshold} | GT == '0')" --fasta-ref ${reference} --missing 'N' --absent 'N' ${vcf} | \
  sed "s/>NC_000962.3 Mycobacterium tuberculosis H37Rv, complete genome/>${sample_id}/g" > ${sample_id}_${params.variant_caller}_mask.fa
 
  # Consensus with no masking of bed_file, still includes filtering by depth/qual thresholds; exclude indels. 
  bcftools consensus --include "(TYPE!='indel' & INFO/DP >= ${params.depth_threshold}) & (QUAL >= ${params.qual_threshold} | GT == '0')" --fasta-ref ${reference} --missing 'N' --absent 'N' ${vcf} | \
  sed "s/>NC_000962.3 Mycobacterium tuberculosis H37Rv, complete genome/>${sample_id}/g" > ${sample_id}_${params.variant_caller}.fa
  """

}
