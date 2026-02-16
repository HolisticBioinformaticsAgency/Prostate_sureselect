process VARDICT_TO_VCF {
  tag "${subject}:${id} (${mode})"
  container 'quay.io/biocontainers/vardict-java:1.8.3--hdfd78af_0'

  publishDir { "${params.outdir_abs}/${subject}/raw_variants" }, mode: 'copy'

  input:
  // (subject, mode, t_sample, n_sample, sample_id, tsv, id)
  tuple val(subject), val(mode), val(t_sample), val(n_sample), val(sample_id), path(tsv), val(id)

  output:
  tuple val(subject), val(id), path("${id}_vardict_raw.vcf.gz"), emit: vcf

  when:
  mode in ['single','paired']

  script:
  // Defensive -N (avoid "null|null")
  def tn = (t_sample && n_sample) ? "${t_sample}|${n_sample}" : "${subject}_T|${subject}_N"

  """
  set -euo pipefail
  echo "[VarDict TO_VCF] Mode: ${mode}" >&2

  if [ "${mode}" = "paired" ]; then
    cat ${tsv} \
      | testsomatic.R \
      | var2vcf_paired.pl -N "${tn}" -f ${params.min_af} \
      > ${id}_vardict_raw.vcf
  else
    cat ${tsv} \
      | teststrandbias.R \
      | var2vcf_valid.pl -N ${sample_id} -E -f ${params.min_af} \
      > ${id}_vardict_raw.vcf
  fi

  gzip -f ${id}_vardict_raw.vcf
  """
}
