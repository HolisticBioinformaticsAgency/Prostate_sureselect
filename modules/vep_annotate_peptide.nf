process VEP_ANNOTATE_PEPTIDE {
  tag "${subject}:${id}"
  container 'quay.io/biocontainers/ensembl-vep:115--pl5321h2a3209d_0'

  // publish per case
  publishDir { "${pub_base}/vep_peptide" }, mode: 'copy'

  input:
  // (pub_base, subject, id, vcf, ref_fa)
  tuple val(pub_base), val(subject), val(id), path(vcf), path(ref_fa)

  output:
  tuple val(pub_base), val(subject), val(id), path("${id}_VEP_${params.vep_species}_${params.vep_assembly}_v${params.vep_version}.vcf.gz"), emit: vepvcf
  tuple val(pub_base), val(subject), val(id), path("${id}_VEP_${params.vep_species}_${params.vep_assembly}_v${params.vep_version}_stats.txt"), emit: stats
  tuple val(pub_base), val(subject), val(id), path("${id}.ref.fa"), path("${id}.mut.fa"), emit: fa
  
  script:
  """
  set -euo pipefail

  OUT_VCF="${id}_VEP_${params.vep_species}_${params.vep_assembly}_v${params.vep_version}.vcf"
  OUT_STATS="${id}_VEP_${params.vep_species}_${params.vep_assembly}_v${params.vep_version}_stats.txt"

  vep \\
    --offline --cache \\
    --species ${params.vep_species} \\
    --assembly ${params.vep_assembly} \\
    --dir_cache ${params.vep_cache_dir} \\
    --fasta ${ref_fa} \\
    --vcf --force_overwrite --format vcf \\
    --buffer_size 50 \\
    --symbol --numbers --mane --hgvs --total_length \\
    --canonical --biotype --protein --ccds --domains --uniprot \\
    --af --af_gnomad --max_af --variant_class \\
    --pick \\
    --plugin ProteinSeqs,${id}.ref.fa,${id}.mut.fa" \\
    --input_file ${vcf} \\
    --output_file \${OUT_VCF} \\
    --stats_text --stats_file \${OUT_STATS}

  gzip -f \${OUT_VCF}
  """
}
