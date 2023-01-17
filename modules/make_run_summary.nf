process SummarizeRun {

  // Parse logs from TrimGalore, Kraken, BWA/Bowtie2
  
  label 'slurm'

  publishDir "${projectDir}/results", mode: "copy", pattern: "pipeline_run_summary.tsv"

  input:
  path reads_list
  path trimming_reports
  path kraken_reports
  path mapping_reports
  path coverage_stats
  path dup_metrics
  path tbprofiler_reports

  output:
  path "pipeline_run_summary.tsv"

  """
  # WORK IN PROGRESS
  touch pipeline_run_summary.tsv
  """

}