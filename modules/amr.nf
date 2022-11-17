process RunAMR {

  // Run AMR prediction tool. Removing for now
  
  label 'slurm'

  publishDir "${projectDir}/results/${batch}/${sample_id}/stats", mode: "copy", pattern: "*_amr.csv"

  input:
  tuple val(sample_id), path(bam), val(batch), val(run)

  output:
  path "${sample_id}_amr.csv"

  """
  mykrobe predict ${bam} tb \
  --output .  \
  --format csv \
  --ploidy haploid \
  --seq ${bam} \
  --panel 201901  # tb 201901 AMR panel based on first line drugs from NEJM-2018 variants (DOI 10.1056/NEJMoa1800474), and second line drugs from Walker 2015 panel.  NC_000962.3
  """

}