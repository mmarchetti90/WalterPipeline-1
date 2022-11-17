process TrimFastQ {

  // Trim reads for quality

  label 'slurm'

  //publishDir "${projectDir}/results/${batch}/${sample_id}/trim", mode: "copy", pattern: "*_val_{1,2}.fq.gz"
  publishDir "${projectDir}/results/${batch}/${sample_id}/trim", mode: "copy", pattern: "*_fastqc.{html,zip}"
  publishDir "${projectDir}/results/${batch}/${sample_id}/trim", mode: "copy", pattern: "*_trimming_report.txt"

  input:
  tuple val(sample_id), path(read1), path(read2), val(batch), val(run)

  output:
  path "*_fastqc.{html,zip}"
  path "*_trimming_report.txt"
  tuple val("${sample_id}"), path("${sample_id}_val_1.fq.gz"), path("${sample_id}_val_2.fq.gz"), val("${batch}"), val("${run}") emit: trimmed_fastq_files

  """
  trim_galore \
  --cores 4 \
  --output_dir . \
  --basename ${sample_id} \
  --nextseq 20
  --fastqc \
  --gzip \
  --paired \
  ${read1} ${read2}
  """

}