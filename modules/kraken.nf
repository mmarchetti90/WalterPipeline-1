process Kraken {

  // Filter reads taxonomically with Kraken
  
  label 'slurm'

  //publishDir "${projectDir}/results/${batch}/${sample_id}/kraken", mode: "copy", pattern: "*_kr_{1,2}.fq.gz"
  //publishDir "${projectDir}/results/${batch}/${sample_id}/kraken", mode: "copy", pattern: "*_kraken.report"
  //publishDir "${projectDir}/results/${batch}/${sample_id}/kraken", mode: "copy", pattern: "*_kraken_stats.csv"

  input:
  path(kraken_db)
  tuple val(sample_id), path(read1), path(read2), val(batch), val(run)

  output:
  //path "*_kraken.report"
  //path "*_kraken_stats.csv"
  tuple val("${sample_id}"), path("${sample_id}_kr_1.fq.gz"), path("${sample_id}_kr_2.fq.gz"), val("${batch}"), val("${run}"), emit: kraken_filtered_files

  """
  # run kraken to taxonomically classify paired-end reads and write output file.
  kraken2 --db . --paired --gzip-compressed --threads \$SLURM_CPUS_ON_NODE --report ${sample_id}_kraken.report --use-names ${read1} ${read2} --output ${sample_id}.out

  # select for reads directly assigned to Mycobacterium genus (G) (taxid 1763), reads assigned directly to Mycobacterium tuberculosis complex (G1) (taxid 77643), and reads assigned to Mycobacterium tuberculosis (S) and children. *this includes reads assigned to Mtb and those assigned to the genus, but not a different species.
  grep -E 'Mycobacterium (taxid 1763)|Mycobacterium tuberculosis' ${sample_id}.out | awk '{print \$2}' > ${sample_id}_reads.list

  # use bbmap to select reads corresponding to taxa of interest.
  filterbyname.sh int=false in1=${read1} in2=${read2} out1=${sample_id}_kr_1.fq.gz out2=${sample_id}_kr_2.fq.gz names=${sample_id}_reads.list include=true overwrite=true
  
  # Summarize Kraken statistics. 
  #${projectDir}/scripts/kraken_stats.sh ${sample_id}_kraken.report
  """

}