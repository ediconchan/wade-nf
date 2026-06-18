//
// GAS (Group A Streptococcus) ALL routine  (SPEC.md section 3)
//   AMR + labware, MLST, EMM, rRNA16S + labware, TOXINS + labware
//
include { MASTERBLASTR as MASTERBLASTR_AMR     } from '../../modules/local/masterblastr'
include { MASTERBLASTR as MASTERBLASTR_TOXINS  } from '../../modules/local/masterblastr'
include { MASTERBLASTR as MASTERBLASTR_RRNA16S } from '../../modules/local/masterblastr'
include { MLST                                 } from '../../modules/local/mlst'
include { EMM                                  } from '../../modules/local/emm'
include { RRNA16S_FORMAT                       } from '../../modules/local/rrna16s_format'
include { LABWARE_GAS_AMR                      } from '../../modules/local/labware_gas_amr'
include { LABWARE_GAS_TOXINS                   } from '../../modules/local/labware_gas_toxins'
include { MLST_SPLIT                           } from '../../modules/local/mlst_split'

workflow GAS {
    take:
    ch_samples   // [ meta(id,org,variable), contig ]
    wade_data    // path to wade-data base dir

    main:
    ch_versions = Channel.empty()

    ref_amr     = file("${wade_data}/GAS/AMR",     checkIfExists: true)
    ref_toxins  = file("${wade_data}/GAS/TOXINS",  checkIfExists: true)
    ref_rrna16s = file("${wade_data}/GAS/rRNA16S", checkIfExists: true)
    ref_mlst    = file("${wade_data}/GAS/MLST",    checkIfExists: true)
    ref_emm     = file("${wade_data}/GAS/EMM",     checkIfExists: true)

    //
    // Per-sample analyses
    //
    MASTERBLASTR_AMR(ch_samples, ref_amr, 'AMR')
    MASTERBLASTR_TOXINS(ch_samples, ref_toxins, 'TOXINS')
    MASTERBLASTR_RRNA16S(ch_samples, ref_rrna16s, 'rRNA16S')
    MLST(ch_samples, ref_mlst, 'MLST')
    EMM(ch_samples, ref_emm)

    ch_versions = ch_versions.mix(
        MASTERBLASTR_AMR.out.versions.first(),
        MASTERBLASTR_TOXINS.out.versions.first(),
        MASTERBLASTR_RRNA16S.out.versions.first(),
        MLST.out.versions.first(),
        EMM.out.versions.first()
    )

    //
    // Collate per-sample rows into one CSV per analysis (replaces R rbind)
    //
    amr_collated = MASTERBLASTR_AMR.out.profile.map { meta, f -> f }
        .collectFile(name: 'output_profile_GAS_AMR.csv', keepHeader: true, sort: { it.text })
    toxins_collated = MASTERBLASTR_TOXINS.out.profile.map { meta, f -> f }
        .collectFile(name: 'output_profile_GAS_TOXINS.csv', keepHeader: true, sort: { it.text })
    rrna16s_collated = MASTERBLASTR_RRNA16S.out.profile.map { meta, f -> f }
        .collectFile(name: 'output_profile_GAS_rRNA16S.csv', keepHeader: true, sort: { it.text })
    mlst_collated = MLST.out.profile.map { meta, f -> f }
        .collectFile(name: 'output_profile_GAS_MLST.csv', keepHeader: true, sort: { it.text })
    emm_collated = EMM.out.profile.map { meta, f -> f }
        .collectFile(name: 'LabWareUpload_GAS_emm.csv', keepHeader: true, sort: { it.text }, storeDir: "${params.outdir}/GAS")

    // Collate extracted-sequence FASTAs across the MASTERBLASTR analyses into one
    // file each (matches the original, which appends all loci to a single
    // output_dna.fasta / output_aa.fasta). Done here to avoid same-named per-task
    // outputs racing in publishDir. sort by content keeps re-runs identical
    // (sort:true would order by random temp-file path -> nondeterministic).
    dna_fasta = MASTERBLASTR_AMR.out.dna.mix(MASTERBLASTR_TOXINS.out.dna, MASTERBLASTR_RRNA16S.out.dna)
        .map { meta, f -> f }
        .collectFile(name: 'output_dna.fasta', sort: { it.text }, storeDir: "${params.outdir}/GAS")
    aa_fasta = MASTERBLASTR_AMR.out.aa.mix(MASTERBLASTR_TOXINS.out.aa, MASTERBLASTR_RRNA16S.out.aa)
        .map { meta, f -> f }
        .collectFile(name: 'output_aa.fasta', sort: { it.text }, storeDir: "${params.outdir}/GAS")

    //
    // LabWare formatters (post-collate, SPEC.md section 4i)
    //
    LABWARE_GAS_AMR(amr_collated)
    LABWARE_GAS_TOXINS(toxins_collated)
    RRNA16S_FORMAT(rrna16s_collated.map { f -> tuple('GAS', f) })
    MLST_SPLIT(mlst_collated.map { f -> tuple('GAS', 'MLST', f) })

    ch_versions = ch_versions.mix(
        LABWARE_GAS_AMR.out.versions,
        LABWARE_GAS_TOXINS.out.versions,
        RRNA16S_FORMAT.out.versions,
        MLST_SPLIT.out.versions
    )

    emit:
    amr      = LABWARE_GAS_AMR.out.labware
    toxins   = LABWARE_GAS_TOXINS.out.labware
    rrna16s  = RRNA16S_FORMAT.out.labware
    mlst     = MLST_SPLIT.out.labware
    emm      = emm_collated
    dna      = dna_fasta
    aa       = aa_fasta
    versions = ch_versions
}
