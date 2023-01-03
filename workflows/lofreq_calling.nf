
// ----------------Workflow---------------- //

include { VariantsLoFreq } from '../modules/variants_lofreq.nf'
include { ConvertVCF } from '../modules/vcf2fasta.nf'
include { AnnotateVCF } from '../modules/annotate_vcf.nf'

workflow LOFREQ {

  take:
  bam_files
	
  main:
  // GATK VARIANT CALLER ------------------ //

  // Channel for genome reference fasta
  reference_fasta = Channel.fromPath(params.reference_fasta_path)

  // Channel for genome reference fasta index
  reference_fasta_index = Channel.fromPath(params.reference_fasta_index_path)

  // Variant calling
  VariantsLoFreq(reference_fasta, reference_fasta_index, bam_files)

  // CONVERTING VCF TO FASTA -------------- //

  ConvertVCF(reference_fasta, VariantsLoFreq.out.lofreq_vcf)

  // ANNOTATE GATK VCF -------------------- //

  AnnotateVCF(reference_fasta, VariantsLoFreq.out.lofreq_vcf)

}