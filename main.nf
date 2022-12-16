#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
INSERT PIPELINE DESCRIPTION
*/

// ----------------Workflow---------------- //

include { TrimFastQ } from './modules/trimgalore.nf'
include { Kraken } from './modules/kraken.nf'
include { QuantTB } from './modules/quanttb.nf'
include { MapReads_BWA } from './modules/map_reads_bwa.nf'
include { MapReads_Bowtie } from './modules/map_reads_bowtie.nf'
include { TbProfiler } from './modules/tb_profiler.nf'
include { RunAMR } from './modules/amr.nf'
include { GATK } from './workflows/gatk_calling.nf'
include { LOFREQ } from './workflows/lofreq_calling.nf'

workflow {

  // CREATING RAW-READS CHANNEL ----------- //

  Channel
    .fromPath(params.reads_list)
    .splitCsv(header: true, sep: '\t')
    .map{row -> tuple(row.sample, row.fastq_1, row.fastq_2, row.batch, row.run)}
    .set{raw_reads}

  // TRIMGALORE --------------------------- //

  TrimFastQ(raw_reads)

  // KRAKEN ------------------------------- //

  // Channel for Kraken2 database
  Channel.fromPath("${params.kraken_database_path}/*{kmer_distrib,k2d,txt,map}")
  .collect()
  .set{kraken_database}

  // Running Kraken2
  Kraken(kraken_database, TrimFastQ.out.trimmed_fastq_files)

  // QUANTTB ------------------------------ //

  //QuantTB(Kraken.out.kraken_filtered_files)

  // MAPPING READS ------------------------ //

  if (params.mapper == "bwa") {

    // MAPPING READS WITH BWA --------------- //

    // Channel for genome reference fasta
    reference_fasta = Channel.fromPath(params.reference_fasta_path)

    // Channel for BWA index
    Channel.fromPath("${params.bwa_index_path}/*{amb,ann,bwt,pac,sa}")
    .collect()
    .set{bwa_index}

    // Mapping and removing duplicates
    MapReads_BWA(reference_fasta, bwa_index, Kraken.out.kraken_filtered_files)

    bam_files = MapReads_BWA.out.bam_files

  }
  else {

    // MAPPING READS WITH BOWTIE2 ----------- //

    // Channel for Bowtie2 index
    Channel.fromPath("${params.bowtie_index_path}/*bt2")
    .collect()
    .set{bowtie_index}

    // Mapping and removing duplicates
    MapReads_Bowtie(bowtie_index, Kraken.out.kraken_filtered_files)

    bam_files = MapReads_Bowtie.out.bam_files

  }

  // TB PROFILER -------------------------- //

  TbProfiler(bam_files)

  // AMR ---------------------------------- //

  RunAMR(bam_files)

  // VARIANT CALLING ---------------------- //

  if (params.variant_caller == "gatk") {

    GATK(bam_files)

  }
  else {

    LOFREQ(bam_files)

  }

}