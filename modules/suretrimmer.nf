process SURETRIMMER {
  tag "${subject}_${sample_id}"
  container ''

  publishDir { "${params.outdir_abs}/${subject}/trimmed" }, mode: 'copy'

  input:
    tuple val(subject), val(sample_id), val(status), path(R1), path(R2)

  output:
  tuple val(subject), val(sample_id), val(status),
    path("${sample_id}.trimmed_R1.fastq.gz"),
    path("${sample_id}.trimmed_R2.fastq.gz")
        

  script:
  """
  bash /projects/vh83/local_software/agent3.0/agent.sh trim -fq1 ${R1} -fq2 ${R2} -v2 -out ./${sample_id}.trimmed
  """
}