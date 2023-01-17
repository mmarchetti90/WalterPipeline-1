process VariantsGATK {

  // Variant calling with GATK
  
  label 'slurm'

  publishDir "${projectDir}/results/${batch}/${sample_id}/vars", mode: "copy", pattern: "*_{gatk.g,gatk_unfilt,gatk_filt}.vcf.gz"

  input:
  each path(reference)
  each path(reference_index)
  each path(dictionary)
  each path(bed)
  each path(bed_index)
  tuple val(sample_id), path(bam), val(batch), val(run)

  output:
  tuple val(sample_id), path("${sample_id}_gatk.g.vcf.gz"), val(batch), val(run), emit: gatk_gvcf
  tuple val(sample_id), path("${sample_id}_gatk_unfilt.vcf.gz"), val(batch), val(run), emit: gatk_vcf_unfilt
  tuple val(sample_id), path("${sample_id}_gatk_filt.vcf.gz"), val(batch), val(run), emit: gatk_vcf_filt

  """
  # Indexing bam
  samtools index ${bam}

  # Call variants with GATK, output GVCF
  # ERC: Reference model emitted with condensed non-variant blocks, i.e. the GVCF format
  gatk --java-options "-Xmx4g" HaplotypeCaller \
  -R ${reference} \
  -ploidy ${params.ploidy} \
  -I ${bam} \
  -ERC GVCF \
  --output-mode EMIT_ALL_CONFIDENT_SITES \
  -O ${sample_id}_gatk.g.vcf.gz

  # GVCF to VCF. Min base quality score is 10 by default. Including non-variant sites in order to differentiate between consensus call and no-call sites.
  gatk --java-options '-Xmx100g' GenotypeGVCFs \
  -R ${reference} \
  --variant ${sample_id}_gatk.g.vcf.gz \
  -ploidy ${params.ploidy} \
  --include-non-variant-sites true \
  --output ${sample_id}_gatk_unfilt.vcf.gz

  # Filter variants by depth and quality score (annotated as PASS or with filter name, not removed).
  gatk VariantFiltration \
  -R ${reference} \
  -V ${sample_id}_gatk_unfilt.vcf.gz \
  -O ${sample_id}_gatk_filt.vcf.gz \
  --filter-name "lowDepth" \
  --filter-expression "DP < ${params.depth_threshold}" \
  --filter-name "lowQual" \
  --filter-expression "QUAL < ${params.qual_threshold}" \
  --mask-name "PPE" \
  --mask ${bed}
  """

}