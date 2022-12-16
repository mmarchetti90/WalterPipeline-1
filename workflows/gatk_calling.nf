
// ----------------Workflow---------------- //

include { VariantsGATK } from '../modules/variants_gatk.nf'
include { ConvertVCF } from '../modules/vcf2fasta.nf'
include { AnnotateVCF } from '../modules/annotate_vcf.nf'

params.caller = "gatk"

workflow GATK {

  take:
  bam_files
	
  main:
	
  // GATK VARIANT CALLER ------------------ //

  // Channel for genome reference fasta
  reference_fasta = Channel.fromPath(params.reference_fasta_path)

  // Channel for genome reference fasta index
  reference_fasta_index = Channel.fromPath(params.reference_fasta_index_path)

  // Channel for GATK dictionary
  gatk_dictionary = Channel.fromPath(params.gatk_dictionary_path)

  // Variant calling
  VariantsGATK(reference_fasta, reference_fasta_index, gatk_dictionary, bam_files)

  // CONVERTING VCF TO FASTA -------------- //

  ConvertVCF(reference_fasta, VariantsGATK.out.gatk_vcf)

  // ANNOTATE GATK VCF -------------------- //

  AnnotateVCF(reference_fasta, VariantsGATK.out.gatk_vcf)

}