nextflow.enable.dsl=2

// ------------------------------------------------------------
// Absolute outdir (stable for publishDir closures)
params.outdir_abs = params.outdir?.startsWith('/') \
  ? params.outdir \
  : "${projectDir}/${params.outdir ?: 'results'}"

// -------------------- Includes --------------------
include { INIT_PARAMS }                    from './modules/init_params.nf'
include { FASTQC }                         from './modules/fastqc.nf'

include { HSMETRICS }                      from './modules/hsmetrics.nf'

include { VARDICT_SINGLE }                 from './modules/vardict.nf'
include { VARDICT_PAIRED }                 from './modules/vardict.nf'

include { SURETRIMMER }                    from './modules/suretrimmer.nf'

include {ALIGN_AND_SORT}                   from './modules/align_and_sort.nf'

include { SET_MATE_INFO; GROUP_READS; GENERATE_CONSENSUS; FGBIO_STATS; MAP_CONSENSUS; INDEX } from './modules/fgbio.nf'
include { VEP_ANNOTATE }                   from './modules/vep_annotate.nf'
include { SNPEFF_ANNOTATE }                from './modules/snpeff_annotate.nf'
include { CLINVAR_ANNOTATE }               from './modules/clinvar_annotate.nf'
include { SOMATIC_FILTER }                 from './modules/somatic_filter.nf'
include { MULTIQC }                        from './modules/multiqc.nf'



// -------------------- Workflow --------------------
workflow {
    // -------------------- Samplesheet → Channels --------------------
def REQUIRED = ['sample','subject','status','fastq_1','fastq_2']

Channel
  .fromPath(params.samplesheet)
  .splitCsv(header:true)
  .map { row ->
    def missing = REQUIRED.findAll { k ->
      !row.containsKey(k) || row[k] == null || row[k].toString().trim() == ''
    }
    if( missing ) {
      throw new IllegalArgumentException("Samplesheet row has missing/empty field(s): ${missing} → ${row}")
    }
    tuple(
      row.sample as String,
      row.subject as String,
      (row.status as String).toLowerCase(),   // 'tumor' | 'normal'
      file(row.fastq_1),
      file(row.fastq_2),
      ((row.sex ?: 'NA') as String)
    )
  }
  .set { ch_sheet }  // (sample, subject, status, R1, R2, sex)

// Reads for QC/BWA: (subject, sample, [R1,R2])
ch_reads = ch_sheet.map { sample, subject, status, r1, r2, sex ->
  tuple(subject, sample, [ r1, r2 ])
}

ch_trim = ch_sheet.map { sample, subject, status, r1, r2, sex ->
  tuple(subject, sample, r1, r2)}

// Sample metadata for joins: (sample, subject, status, sex)
ch_meta  = ch_sheet.map { sample, subject, status, r1, r2, sex ->
  tuple(sample, subject, status, sex)
}

// Reference & BED singletons
Channel.of( file(params.reference)     ).set { ch_ref_fa }
Channel.of( file(params.reference_fai) ).set { ch_ref_fai }
Channel.of( file(params.bed)           ).set { ch_bed }

// Also pass original absolute FASTA path (for BWA sidecars)
Channel
  .of( params.reference )
  .map { new File(it as String).getAbsolutePath() }
  .set { ch_ref_src_abs }
  // ---- Param check ---------
  INIT_PARAMS()

  // ---------- QC ----------
  FASTQC( ch_reads )

  // ---------- trim reads ------------
  ch_trim_out = SURETRIMMER ( ch_trim )
  ch_align_in_alpha = ch_trim_out.map { sample, subject, r1, r2 -> 
    tuple(sample, subject, [r1, r2])}

  // ---------- Align + Sort ----------
  ch_align_in   = ch_align_in_alpha.combine(ch_ref_fa).combine(ch_ref_src_abs)
  ch_bam_sorted = ALIGN_AND_SORT( ch_align_in )

  // ---------- FGBIO ----------
  
  ch_fgbio_out1 = SET_MATE_INFO(ch_bam_sorted, params.reference)
  ch_fgbio_out2 = GROUP_READS (ch_fgbio_out1)
  ch_fgbio_out3 = GENERATE_CONSENSUS (ch_fgbio_out2)
  ch_fgbio_out4 = FGBIO_STATS (ch_fgbio_out3)
  ch_fgbio_out5 = MAP_CONSENSUS ( ch_fgbio_out3.combine(ch_ref_fa).combine(ch_ref_src_abs))
  ch_bam = INDEX ( ch_fgbio_out5 )
  
  ch_fgbio_out5.view()
  
  
 
  // ---------- Attach metadata ----------
  ch_bam_meta = ch_bam.collect()
    .map  { sub, sample, bam, bai -> tuple(sample, tuple(sub, bam, bai)) }
    .join ( ch_meta )
    .map  { sample, left, subject, status, sex ->
      tuple(subject, sample, status, sex, left[1], left[2])
      
    }


  // ---------- HsMetrics ----------
  def ch_hs_in = ch_bam.combine(ch_bed).combine(ch_ref_fa).combine(ch_ref_fai)
    .map { sub, sample, bam, bai, bed, ref_fa, ref_fai -> tuple(sub, sample, bam, bai, bed, ref_fa, ref_fai) }
  HSMETRICS( ch_hs_in )

 // ---------- Build per-subject groups (deterministic) ----------
  def by_subject = ch_bam_meta
    .map { sub, sample, status, sex, bam, bai -> tuple(sub, tuple(sample, status, bam, bai)) }
    .groupTuple()
    .map { sub, recs -> tuple(sub, recs.sort { a, b -> a[0] <=> b[0] }) } // stable order
    

  // ---------- Compute cases and publishing base in one pass ----------
  def ch_cases_pub = by_subject.collect().flatMap { sub, recs ->
    // rec: [sample, status, bam, bai]
    def tumors  = recs.findAll { it[1] == 'tumor'  }.collect { [ it[0], it[2], it[3], 'tumor'  ] } // [id,bam,bai,status]
    def normals = recs.findAll { it[1] == 'normal' }.collect { [ it[0], it[2], it[3], 'normal' ] }

    // number of cases for this subject
    def ncases = (tumors && normals) ? (tumors.size() * normals.size()) : recs.size()

    def out = []
    if( tumors && normals ) {
      tumors.each { t ->
        normals.each { n ->
          def case_id = "${t[0]}_${n[0]}"
          def base    = (ncases > 1) ? "${params.outdir_abs}/${sub}/${case_id}"
                                    : "${params.outdir_abs}/${sub}"
          out << tuple(sub, case_id, 'paired', [t, n], base)
        }
      }
    } else {
      recs.each { r ->
        def case_id = r[0]
        def base    = (ncases > 1) ? "${params.outdir_abs}/${sub}/${case_id}"
                                  : "${params.outdir_abs}/${sub}"
        out << tuple(sub, case_id, 'single', [ [ r[0], r[2], r[3], r[1] ] ], base)
      }
    }
    out
  }
  
  

  // ---------- VarDict ----------
  def ch_vardict_single_in = ch_cases_pub
    .filter { sub, case_id, mode, samples, pub_base -> mode == 'single' }
    .combine(ch_ref_fa).combine(ch_ref_fai).combine(ch_bed)
    .map { sub, case_id, mode, samples, pub_base, ref_fa, ref_fai, bed ->
      def s = samples[0]
      tuple(pub_base, sub, s[0], s[1], s[2], ref_fa, ref_fai, bed)
    }

  def ch_vardict_paired_in = ch_cases_pub
    .filter { sub, case_id, mode, samples, pub_base -> mode == 'paired' }
    .combine(ch_ref_fa).combine(ch_ref_fai).combine(ch_bed)
    .flatMap { sub, case_id, mode, samples, pub_base, ref_fa, ref_fai, bed ->
      def tumors  = samples.findAll { it[3] == 'tumor'  }
      def normals = samples.findAll { it[3] == 'normal' }
      tumors.collectMany { t ->
        normals.collect { n ->
          tuple(pub_base, sub, case_id, t[0], t[1], t[2], n[0], n[1], n[2], ref_fa, ref_fai, bed)
        }
      }
    }

  def ch_vardict_single_vcf = VARDICT_SINGLE( ch_vardict_single_in )
    .map { pub_base, sub, id, vcf -> tuple(pub_base, sub, id, 'single', vcf) }

  def ch_vardict_paired_vcf = VARDICT_PAIRED( ch_vardict_paired_in )
    .map { pub_base, sub, case_id, vcf -> tuple(pub_base, sub, case_id, 'paired', vcf) }

  def ch_variants_vcf = ch_vardict_single_vcf.mix( ch_vardict_paired_vcf )


  // ---------- Annotation ----------
  def ch_vep_in = ch_variants_vcf.combine(ch_ref_fa)
    .map { pub_base, sub, id, mode, vcf, ref_fa ->
      tuple(pub_base, sub, id, file(vcf), file(ref_fa))
    }
  def (ch_vep_vcf, ch_vep_stats) = VEP_ANNOTATE( ch_vep_in )
/*
  def ch_mane_dir    = Channel.of( file(params.mane_dir) )
  def ch_snpeff_core = SNPEFF_ANNOTATE(
    ch_variants_vcf.combine(ch_mane_dir).map { pub_base, sub, id, mode, vcf, mane ->
      tuple(pub_base, sub, id, mode, vcf, mane)
    }
  )
  def ch_clinvar_vcf = Channel.of( file(params.clinvar_vcf) )
  def ch_clinvar_out = CLINVAR_ANNOTATE(
    ch_snpeff_core.combine(ch_clinvar_vcf).map { pub_base, sub, id, mode, core_vcf, clin ->
      tuple(pub_base, sub, id, mode, core_vcf, clin)
    }
  )
*/
/*
  def ch_annot_for_all = (params.pytmb_annot == 'snpeff')
    ? ch_clinvar_out.map { pub_base, sub, id, mode, vcf -> tuple(pub_base, sub, id, mode, vcf) }
    : ch_vep_vcf.combine( ch_variants_vcf.map{ pb, s, i, m, v -> tuple(s,i,m) } )
                 .map { vep_pb_sub_id_vcf, sub_id_mode ->
                   def (pb, s1, i1, vepvcf) = vep_pb_sub_id_vcf
                   def (s2, i2, mode)      = sub_id_mode
                   assert s1==s2 && i1==i2
                   tuple(pb, s1, i1, mode, vepvcf)
                 }

  // Stop normal-only after annotation
  def ch_cases_has_tumor = ch_cases_pub.map { sub, case_id, mode, samples, pub_base ->
    tuple("${sub}::${case_id}", sub, case_id, samples.any{ it[3]=='tumor' })
  }
  def ch_annot_keyed    = ch_annot_for_all.map { pub, sub, id, mode, vcf -> tuple("${sub}::${id}", pub, sub, id, mode, vcf) }
  def ch_anno_with_flag = ch_annot_keyed.join( ch_cases_has_tumor )
    .map { key, pub, sub, id, mode, vcf, sub2, case_id, hasTumor ->
      tuple(pub, sub, id, mode, vcf, hasTumor)
    }

  def ch_for_tumor = ch_anno_with_flag
    .filter { pub, sub, id, mode, vcf, hasTumor -> hasTumor }
    .map    { pub, sub, id, mode, vcf, hasTumor -> tuple(pub, sub, id, mode, vcf) }

  def ch_somatic_filtered_vcf = SOMATIC_FILTER( ch_for_tumor )

  def ch_vep_peptide_in = ch_somatic_filtered_vcf.combine(ch_ref_fa)
    .map { pub, sub, id, vcf, ref -> tuple(pub, sub, id, vcf, ref)
    }

  VEP_ANNOTATE_PEPTIDE(ch_vep_peptide_in)

  def ch_cases_keyed = ch_cases_pub.map { sub, case_id, mode, samples, pub_base ->
    tuple("${sub}::${case_id}", sub, case_id, mode, samples, pub_base)
  }
  def ch_somatic_key = ch_somatic_filtered_vcf.map { pub, sub, id, vcf ->
    tuple("${sub}::${id}", pub, sub, id, vcf)
  }

  
  

  // ---------- MultiQC ----------
  def ch_hsmetrics_files = HSMETRICS.out.hs.map { sub, sid, hs -> hs }
  def ch_multiqc_inputs  = ch_dedup_metrics.mix(ch_hsmetrics_files).collect()
  MULTIQC( ch_multiqc_inputs )
 

workflow.onComplete {
    def outdir = params.outdir_abs
    def cmd = """
      set -euo pipefail
      shopt -s nullglob

      for d in "${outdir}"/* ; do
        [ -d "\$d" ] || continue
        base="\$(basename "\$d")"
        # Skip non-subject folders at top level
        case "\$base" in
          reference|multiqc) continue ;;
        esac

        dest="\${d%/}/vcf"
        mkdir -p "\$dest"

        # Gather all VCFs & indexes anywhere under the subject dir,
        # but don't recurse into the destination vcf/ we're filling.
        find "\$d" -type f \\
          \\( -name "*.vcf" -o -name "*.vcf.gz" -o -name "*.vcf.tbi" -o -name "*.vcf.gz.tbi" -o -name "*.csi" \\) \\
          -not -path "\$dest/*" \\
          -exec cp -f {} "\$dest/" \\;
      done
    """
    def p = ["bash","-lc", cmd].execute()
    p.consumeProcessOutput(System.out, System.err)
    def rc = p.waitFor()
    if( rc != 0 ) log.warn "VCF gather post-step exited with code ${rc}"
}
*/

}

