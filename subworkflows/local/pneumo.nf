//
// PNEUMO (Streptococcus pneumoniae) ALL routine  (SPEC.md section 3)
//   AMR_ALL (AMR + rRNA23S) + labware, MLST, VIRULENCE + labware, SERO
//
include { MASTERBLASTR as MASTERBLASTR_AMR       } from '../../modules/local/masterblastr'
include { MASTERBLASTR as MASTERBLASTR_VIRULENCE } from '../../modules/local/masterblastr'
include { MLST                                   } from '../../modules/local/mlst'
include { SEROTYPE                               } from '../../modules/local/serotype'
include { RRNA23S                                } from '../../modules/local/rrna23s'
include { LABWARE_PNEUMO_AMR                     } from '../../modules/local/labware_pneumo_amr'
include { LABWARE_PNEUMO_VIRULENCE               } from '../../modules/local/labware_pneumo_virulence'
include { MLST_SPLIT                             } from '../../modules/local/mlst_split'

workflow PNEUMO {
    take:
    ch_samples   // [ meta(id,org,variable), contig ]
    ch_vcf       // [ meta, vcf ]   (vcf may be [])
    wade_data

    main:
    ch_versions = Channel.empty()

    ref_amr       = file("${wade_data}/PNEUMO/AMR",       checkIfExists: true)
    ref_virulence = file("${wade_data}/PNEUMO/VIRULENCE", checkIfExists: true)
    ref_mlst      = file("${wade_data}/PNEUMO/MLST",      checkIfExists: true)
    ref_sero      = file("${wade_data}/PNEUMO/SERO",      checkIfExists: true)

    MASTERBLASTR_AMR(ch_samples, ref_amr, 'AMR')
    MASTERBLASTR_VIRULENCE(ch_samples, ref_virulence, 'VIRULENCE')
    RRNA23S(ch_vcf)
    MLST(ch_samples, ref_mlst, 'MLST')
    SEROTYPE(ch_samples, ref_sero)

    ch_versions = ch_versions.mix(
        MASTERBLASTR_AMR.out.versions.first(),
        MASTERBLASTR_VIRULENCE.out.versions.first(),
        RRNA23S.out.versions.first(),
        MLST.out.versions.first(),
        SEROTYPE.out.versions.first()
    )

    // Collate per-sample rows
    amr_collated = MASTERBLASTR_AMR.out.profile.map { meta, f -> f }
        .collectFile(name: 'output_profile_PNEUMO_AMR.csv', keepHeader: true, sort: { it.text })
    virulence_collated = MASTERBLASTR_VIRULENCE.out.profile.map { meta, f -> f }
        .collectFile(name: 'output_profile_PNEUMO_VIRULENCE.csv', keepHeader: true, sort: { it.text })
    rrna23s_collated = RRNA23S.out.profile.map { meta, f -> f }
        .collectFile(name: 'output_profile_PNEUMO_rRNA23S.csv', keepHeader: true, sort: { it.text })
    mlst_collated = MLST.out.profile.map { meta, f -> f }
        .collectFile(name: 'output_profile_PNEUMO_MLST.csv', keepHeader: true, sort: { it.text })
    sero_collated = SEROTYPE.out.profile.map { meta, f -> f }
        .collectFile(name: 'LabWareUpload_PNEUMO_SEROTYPE.csv', keepHeader: true, sort: { it.text }, storeDir: "${params.outdir}/PNEUMO")

    dna_fasta = MASTERBLASTR_AMR.out.dna.mix(MASTERBLASTR_VIRULENCE.out.dna)
        .map { meta, f -> f }
        .collectFile(name: 'output_dna.fasta', sort: { it.text }, storeDir: "${params.outdir}/PNEUMO")
    aa_fasta = MASTERBLASTR_AMR.out.aa.mix(MASTERBLASTR_VIRULENCE.out.aa)
        .map { meta, f -> f }
        .collectFile(name: 'output_aa.fasta', sort: { it.text }, storeDir: "${params.outdir}/PNEUMO")

    // AMR_ALL chain: AMR profile + 23S counts -> pneumo AMR labware
    LABWARE_PNEUMO_AMR(amr_collated, rrna23s_collated, ref_amr)
    LABWARE_PNEUMO_VIRULENCE(virulence_collated)
    MLST_SPLIT(mlst_collated.map { f -> tuple('PNEUMO', 'MLST', f) })

    ch_versions = ch_versions.mix(
        LABWARE_PNEUMO_AMR.out.versions,
        LABWARE_PNEUMO_VIRULENCE.out.versions,
        MLST_SPLIT.out.versions
    )

    emit:
    amr       = LABWARE_PNEUMO_AMR.out.labware
    virulence = LABWARE_PNEUMO_VIRULENCE.out.labware
    mlst      = MLST_SPLIT.out.labware
    sero      = sero_collated
    versions  = ch_versions
}
