process VARDICT_SINGLE_RAW {
  tag "${subject}_${sample_id}"
  container 'quay.io/biocontainers/vardict-java:1.8.3--hdfd78af_0'
  env = [ 'JAVA_TOOL_OPTIONS': '-Xms2g -Xmx20g' ]

  // per-subject outputs
  publishDir { "${params.outdir_abs}/${subject}/vardict_single" }, mode: 'copy'

  input:
  // from wiring: (subject, sample, bam, bai, ref_fa, ref_fai, bed)
  tuple val(subject), val(sample_id), path(bam), path(bai), path(ref_fa), path(ref_fai), path(bed)

  output:
  // keep subject for downstream joins
  tuple val(subject), val(sample_id), path("${sample_id}_vardict.tsv"), emit: tsv

  script:
  """
  set -euo pipefail

  # Determine BED gene column (3-col vs 4+-col BED)
  GCOL=\$(awk 'BEGIN{c=0} !/^#/ && NF{print (NF>=4?4:0); exit}' ${bed})
  echo "[VarDict SINGLE] Using -g \${GCOL}" >&2

  vardict-java \\
    -G ${ref_fa} \\
    -N ${sample_id} \\
    -b ${bam} \\
    -th ${task.cpus} \\
    -f ${params.min_af} -Q 20 -z -c 1 -S 2 -E 3 -g \${GCOL} ${bed} \\
    > ${sample_id}_vardict.tsv

  test -s ${sample_id}_vardict.tsv || { echo "[VarDict SINGLE] Empty TSV for ${sample_id}" >&2; exit 2; }
  """
}
