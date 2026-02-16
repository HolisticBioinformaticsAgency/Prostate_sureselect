process CLINVAR_ANNOTATE {
  tag "${subject}:${id}"
  container '/fs04/scratch2/cg90/liem/projects/nextflow/molbi_pipeline/containers/snpsift_4.3.sif'

  publishDir { "${pub_base}/snpeff" }, mode: 'copy'

  input:
  // (pub_base, subject, id, mode, mane_vcf, clinvar_vcf)
  tuple val(pub_base), val(subject), val(id), val(mode), path(mane_vcf), path(clinvar_vcf)

  output:
  // (pub_base, subject, id, mode, <id>_MANE_* .vcf.gz)
  tuple val(pub_base), val(subject), val(id), val(mode), path("${id}_MANE_*.vcf.gz"), emit: snpeffvcf

  script:
  """
  set -euo pipefail

  clin_stem=\$(basename "${clinvar_vcf}")
  clin_stem="\${clin_stem%.gz}"
  clin_stem="\${clin_stem%.vcf}"

  OUT_PREFIX="${id}_MANE_\${clin_stem}"
  OUT_VCF="\${OUT_PREFIX}.vcf"
  OUT_GZ="\${OUT_VCF}.gz"

  set +e
  SnpSift annotate -tabix "${clinvar_vcf}" "${mane_vcf}" > "\${OUT_VCF}"
  rc=\$?
  set -e

  if [[ \$rc -ne 0 ]]; then
    echo "[CLINVAR_ANNOTATE] WARN: -tabix annotate failed. Falling back..." >&2
    gunzip -c "${clinvar_vcf}" > "${id}_clinvar_db.vcf"
    SnpSift annotate -sorted "${id}_clinvar_db.vcf" "${mane_vcf}" > "\${OUT_VCF}"
  fi

  gzip -f "\${OUT_VCF}"
  """
}
