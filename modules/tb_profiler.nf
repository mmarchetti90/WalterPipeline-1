process TbProfiler {

  // TB Profiler to assign sub-lineage
  
  label 'slurm'

  publishDir "${projectDir}/results/${batch}/${sample_id}/stats", mode: "copy", pattern: "*_lineageSpo.csv"

  input:
  tuple val(sample_id), path(bam), val(batch), val(run)

  output:
  path "${sample_id}_lineageSpo.csv"

  """
  tb-profiler profile --bam ${bam} --prefix ${sample_id} --dir . --csv --spoligotype
  mv ${sample_id}.results.csv ${sample_id}_lineageSpo.csv
  """

}