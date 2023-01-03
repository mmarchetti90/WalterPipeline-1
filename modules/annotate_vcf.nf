process AnnotateVCF {

  // Annotate VCF for variant examination
  
  label 'slurm'

  publishDir "${projectDir}/results/${batch}/${sample_id}/vars", mode: "copy", pattern: "*_${params.variant_caller}_ann.vcf.gz"

  input:
  each path(reference)
  tuple val(sample_id), path(vcf), val(batch), val(run)

  output:
  path "${sample_id}_${params.variant_caller}_ann.vcf.gz"

  """
  # Rename Chromosome to be consistent with snpEff/Ensembl genomes.
  zcat ${vcf} | sed 's/NC_000962.3/Chromosome/g' | bgzip > ${sample_id}_renamed.vcf.gz
  tabix ${sample_id}_renamed.vcf.gz

  # Run snpEff and then rename Chromosome.
  java -jar -Xmx8g ${params.snpeff} eff ${params.snpeff_db} ${sample_id}_renamed.vcf.gz -dataDir ${params.snpeff_datapath} -noStats -no-downstream -no-upstream -canon | sed 's/Chromosome/NC_000962.3/g' > ${sample_id}_tmp.vcf.gz

  # Also use bed file to annotate vcf, zip.
  bcftools annotate -a ${params.bed_path} -h ${params.vcf_header} -c CHROM,FROM,TO,FORMAT/PPE ${sample_id}_tmp.vcf.gz | bgzip > ${sample_id}_${params.variant_caller}_ann.vcf.gz
  """

}