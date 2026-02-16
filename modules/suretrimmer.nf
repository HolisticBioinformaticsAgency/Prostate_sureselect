process SURETRIMMER {
  tag "${subject}_${sample_id}"
  container ''

  publishDir { "${params.outdir_abs}/${subject}/trimmed" }, mode: 'copy'

  input:
  // (subject, sample_id, [R1,R2], staged FASTA, ABSOLUTE original FASTA path string)
  tuple val(subject), val(sample_id), path(R1), path(R2)

  output:
  tuple val(subject), val(sample_id),
        path("${sample_id}.trimmed_R1.fastq.gz"),
        path("${sample_id}.trimmed_R2.fastq.gz"),
        

  script:
  """
  bash /projects/vh83/local_software/agent3.0/agent.sh trim -fq1 ${R1} -fq2 ${R2} -v2 -out ./${sample_id}.trimmed
  """
}