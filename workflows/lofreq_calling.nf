
// ----------------Workflow---------------- //

include { VariantsLoFreq } from '../modules/variants_lofreq.nf'
include { ConvertVCF } from '../modules/vcf2fasta.nf'
include { AnnotateVCF } from '../modules/annotate_vcf.nf'

params.caller = "lofreq"

workflow LOFREQ {

  take:
  bam_files
	
  main:
	
  // GATK VARIANT CALLER ------------------ //

  // Channel for genome reference fasta
  reference_fasta = Channel.fromPath(params.reference_fasta_path)

  // Variant calling
  VariantsLoFreq(reference_fasta, bam_files)

  // CONVERTING VCF TO FASTA -------------- //

  ConvertVCF(reference_fasta, VariantsGATK.out.lofreq_vcf)

  // ANNOTATE GATK VCF -------------------- //

  AnnotateVCF(reference_fasta, VariantsGATK.out.lofreq_vcf)

}