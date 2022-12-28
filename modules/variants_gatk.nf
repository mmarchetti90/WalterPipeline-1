process VariantsGATK {

  // Variant calling with GATK
  
  label 'slurm'

  publishDir "${projectDir}/results/${batch}/${sample_id}/vars", mode: "copy", pattern: "*_gatk.vcf.gz"

  input:
  each path(reference)
  each path(reference_index)
  each path(dictionary)
  tuple val(sample_id), path(bam), val(batch), val(run)

  output:
  tuple val(sample_id), path("${sample_id}_gatk.vcf.gz"), val(batch), val(run), emit: gatk_vcf

  """
  # Indexing bam
  samtools index ${bam}

  # Call variants with GATK 4.1, output GVCF
  gatk --java-options "-Xmx100g" HaplotypeCaller \
  -R ${reference} \
  -ploidy 1 \
  -I ${bam} \
  -ERC GVCF \
  -O ${sample_id}_gatk.g.vcf.gz

  # Index gvcf 
  gatk IndexFeatureFile \
  -I ${sample_id}_gatk.g.vcf.gz

  # GVCF to VCF. Min base quality score is 10 by default.
  gatk --java-options '-Xmx100g' GenotypeGVCFs \
  -R ${reference} \
  --variant ${sample_id}_gatk.g.vcf.gz \
  -ploidy 1 \
  --include-non-variant-sites true \
  --output ${sample_id}_gatk.vcf.gz
  """

}
