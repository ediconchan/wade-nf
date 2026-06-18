//
// GBS (Group B Streptococcus) ALL routine  (SPEC.md section 3)
//   AMR + labware, MLST, rRNA16S + labware, SERO
//
include { MASTERBLASTR as MASTERBLASTR_AMR     } from '../../modules/local/masterblastr'
include { MASTERBLASTR as MASTERBLASTR_RRNA16S } from '../../modules/local/masterblastr'
include { MLST                                 } from '../../modules/local/mlst'
include { SEROTYPE                             } from '../../modules/local/serotype'
include { RRNA16S_FORMAT                       } from '../../modules/local/rrna16s_format'
include { LABWARE_GBS_AMR                      } from '../../modules/local/labware_gbs_amr'
include { MLST_SPLIT                           } from '../../modules/local/mlst_split'

workflow GBS {
    take:
    ch_samples   // [ meta(id,org,variable), contig ]
    wade_data

    main:
    ch_versions = Channel.empty()

    ref_amr     = file("${wade_data}/GBS/AMR",     checkIfExists: true)
    ref_rrna16s = file("${wade_data}/GBS/rRNA16S", checkIfExists: true)
    ref_mlst    = file("${wade_data}/GBS/MLST",    checkIfExists: true)
    ref_sero    = file("${wade_data}/GBS/SERO",    checkIfExists: true)

    MASTERBLASTR_AMR(ch_samples, ref_amr, 'AMR')
    MASTERBLASTR_RRNA16S(ch_samples, ref_rrna16s, 'rRNA16S')
    MLST(ch_samples, ref_mlst, 'MLST')
    SEROTYPE(ch_samples, ref_sero)

    ch_versions = ch_versions.mix(
        MASTERBLASTR_AMR.out.versions.first(),
        MASTERBLASTR_RRNA16S.out.versions.first(),
        MLST.out.versions.first(),
        SEROTYPE.out.versions.first()
    )

    // Collate per-sample rows
    amr_collated = MASTERBLASTR_AMR.out.profile.map { meta, f -> f }
        .collectFile(name: 'output_profile_GBS_AMR.csv', keepHeader: true, sort: { it.text })
    rrna16s_collated = MASTERBLASTR_RRNA16S.out.profile.map { meta, f -> f }
        .collectFile(name: 'output_profile_GBS_rRNA16S.csv', keepHeader: true, sort: { it.text })
    mlst_collated = MLST.out.profile.map { meta, f -> f }
        .collectFile(name: 'output_profile_GBS_MLST.csv', keepHeader: true, sort: { it.text })
    sero_collated = SEROTYPE.out.profile.map { meta, f -> f }
        .collectFile(name: 'LabWareUpload_GBS_SEROTYPE.csv', keepHeader: true, sort: { it.text }, storeDir: "${params.outdir}/GBS")

    dna_fasta = MASTERBLASTR_AMR.out.dna.mix(MASTERBLASTR_RRNA16S.out.dna)
        .map { meta, f -> f }
        .collectFile(name: 'output_dna.fasta', sort: { it.text }, storeDir: "${params.outdir}/GBS")
    aa_fasta = MASTERBLASTR_AMR.out.aa.mix(MASTERBLASTR_RRNA16S.out.aa)
        .map { meta, f -> f }
        .collectFile(name: 'output_aa.fasta', sort: { it.text }, storeDir: "${params.outdir}/GBS")

    LABWARE_GBS_AMR(amr_collated)
    RRNA16S_FORMAT(rrna16s_collated.map { f -> tuple('GBS', f) })
    MLST_SPLIT(mlst_collated.map { f -> tuple('GBS', 'MLST', f) })

    ch_versions = ch_versions.mix(
        LABWARE_GBS_AMR.out.versions,
        RRNA16S_FORMAT.out.versions,
        MLST_SPLIT.out.versions
    )

    emit:
    amr      = LABWARE_GBS_AMR.out.labware
    rrna16s  = RRNA16S_FORMAT.out.labware
    mlst     = MLST_SPLIT.out.labware
    sero     = sero_collated
    versions = ch_versions
}
