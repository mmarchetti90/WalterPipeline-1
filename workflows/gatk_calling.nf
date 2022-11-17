
// ----------------Workflow---------------- //

include { VariantsGATK } from './modules/variants_gatk.nf'
include { ConvertVCF } from './modules/vcf2fasta.nf'
include { AnnotateVCF } from './modules/annotate_vcf.nf'

workflow GATK {

  take:
  bam_files
	
	// GATK VARIANT CALLER ------------------ //

    // Channel for genome reference fasta
    reference_fasta = Channel.fromPath(params.reference_fasta)

    // Channel for GATK dictionary
    gatk_dictionary = Channel.fromPath(params.gatk_dictionary)

    // Variant calling
    VariantsGATK(reference_fasta, gatk_dictionary, bam_files)

    // CONVERTING VCF TO FASTA -------------- //

    ConvertVCF("gatk", reference_fasta, VariantsGATK.out.gatk_vcf)

    // ANNOTATE GATK VCF -------------------- //

    AnnotateVCF("gatk", VariantsGATK.out.gatk_vcf)

}