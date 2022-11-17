
// ----------------Workflow---------------- //

include { VariantsLoFreq } from './modules/variants_lofreq.nf'
include { ConvertVCF } from './modules/vcf2fasta.nf'
include { AnnotateVCF } from './modules/annotate_vcf.nf'

workflow LOFREQ {

  take:
  bam_files
	
	// GATK VARIANT CALLER ------------------ //

    // Channel for genome reference fasta
    reference_fasta = Channel.fromPath(params.reference_fasta)

    // Variant calling
    VariantsLoFreq(reference_fasta, bam_files)

    // CONVERTING VCF TO FASTA -------------- //

    ConvertVCF("lofreq", reference_fasta, VariantsGATK.out.lofreq_vcf)

    // ANNOTATE GATK VCF -------------------- //

    AnnotateVCF("lofreq", VariantsGATK.out.lofreq_vcf)

}