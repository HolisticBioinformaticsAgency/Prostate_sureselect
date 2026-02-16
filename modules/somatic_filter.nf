process SOMATIC_FILTER {
  tag "${subject}:${id}"
  container '/fs04/scratch2/cg90/liem/projects/nextflow/molbi_pipeline/containers/snpsift_4.3.sif'

  publishDir { "${pub_base}/somatic" }, mode: 'copy'

  input:
  // (pub_base, subject, id, mode, vcf)
  tuple val(pub_base), val(subject), val(id), val(mode), path(vcf)

  output:
  tuple val(pub_base), val(subject), val(id), path("${id}_somatic_filtered.vcf.gz"), emit: somatic

  script:
  """
  set -euo pipefail

  if [ "${mode}" = "paired" ]; then
    echo "[SomaticFilter] Paired mode: filtering STATUS" >&2
    SnpSift filter "(STATUS = 'StrongLOH') | (STATUS = 'StrongSomatic')" "${vcf}" \
      > "${id}_somatic_filtered.vcf"
  else
    echo "[SomaticFilter] Single mode: skipping STATUS filter" >&2
    zcat -f "${vcf}" > "${id}_somatic_filtered.vcf"
  fi

  gzip -f "${id}_somatic_filtered.vcf"
  """
}
