process VARDICT_PAIRED_RAW {
  tag "${subject} (T:${t_sample}|N:${n_sample})"
  container 'quay.io/biocontainers/vardict-java:1.8.3--hdfd78af_0'
  env = [ 'JAVA_TOOL_OPTIONS': '-Xms2g -Xmx20g' ]

  // per-subject outputs
  publishDir { "${params.outdir_abs}/${subject}/vardict_paired" }, mode: 'copy'

  input:
  // from wiring: (subject, t_sample, t_bam, t_bai, n_sample, n_bam, n_bai, ref_fa, ref_fai, bed)
  tuple val(subject),
        val(t_sample), path(t_bam), path(t_bai),
        val(n_sample), path(n_bam), path(n_bai),
        path(ref_fa), path(ref_fai),
        path(bed)

  output:
  // keep subject for downstream joins; use subject as the case_id
  tuple val(subject), path("${subject}_vardict.tsv"), emit: tsv

  script:
  """
  set -euo pipefail

  # Determine BED gene column (3-col vs 4+-col BED)
  GCOL=\$(awk 'BEGIN{c=0} !/^#/ && NF{print (NF>=4?4:0); exit}' ${bed})
  echo "[VarDict PAIRED] Using -g \${GCOL}" >&2

  vardict-java \\
    -G ${ref_fa} \\
    -N ${subject} \\
    -b "${t_bam}|${n_bam}" \\
    -th ${task.cpus} \\
    -f ${params.min_af} -Q 20 -z -c 1 -S 2 -E 3 -g \${GCOL} ${bed} \\
    > ${subject}_vardict.tsv

  test -s ${subject}_vardict.tsv || { echo "[VarDict PAIRED] Empty TSV for ${subject}" >&2; exit 2; }
  """
}
