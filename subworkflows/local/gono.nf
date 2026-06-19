//
// GONO (Neisseria gonorrhoeae) ALL routine  (SPEC.md section 3)
//   AMR_ALL chain (MASTER AMR + NGSTAR_AMR mutations + rRNA23S) + labware,
//   MLST, NGSTAR, NGMAST
//
include { MASTERBLASTR as MASTERBLASTR_AMR } from '../../modules/local/masterblastr'
include { MLST as MLST_NGSTAR_AMR          } from '../../modules/local/mlst'
include { MLST as MLST_MLST                } from '../../modules/local/mlst'
include { MLST as MLST_NGSTAR              } from '../../modules/local/mlst'
include { MLST as MLST_NGMAST              } from '../../modules/local/mlst'
include { RRNA23S                          } from '../../modules/local/rrna23s'
include { LABWARE_GONO_AMR                 } from '../../modules/local/labware_gono_amr'
include { MLST_SPLIT                       } from '../../modules/local/mlst_split'

workflow GONO {
    take:
    ch_samples   // [ meta(id,org,variable), contig ]
    ch_vcf       // [ meta, vcf ]
    wade_data

    main:
    ch_versions = Channel.empty()

    ref_amr    = file("${wade_data}/GONO/AMR",    checkIfExists: true)
    ref_ngstar = file("${wade_data}/GONO/NGSTAR", checkIfExists: true)
    ref_mlst   = file("${wade_data}/GONO/MLST",   checkIfExists: true)
    ref_ngmast = file("${wade_data}/GONO/NGMAST", checkIfExists: true)

    MASTERBLASTR_AMR(ch_samples, ref_amr, 'AMR_ALL')
    MLST_NGSTAR_AMR(ch_samples, ref_ngstar, 'NGSTAR_AMR')
    RRNA23S(ch_vcf)
    MLST_MLST(ch_samples, ref_mlst, 'MLST')
    MLST_NGSTAR(ch_samples, ref_ngstar, 'NGSTAR')
    MLST_NGMAST(ch_samples, ref_ngmast, 'NGMAST')

    ch_versions = ch_versions.mix(
        MASTERBLASTR_AMR.out.versions.first(),
        MLST_NGSTAR_AMR.out.versions.first(),
        RRNA23S.out.versions.first(),
        MLST_MLST.out.versions.first(),
        MLST_NGSTAR.out.versions.first(),
        MLST_NGMAST.out.versions.first()
    )

    // Collate per-sample rows
    amr_collated = MASTERBLASTR_AMR.out.profile.map { meta, f -> f }
        .collectFile(name: 'output_profile_GONO_AMR.csv', keepHeader: true, sort: { it.text })
    ngstar_mut_collated = MLST_NGSTAR_AMR.out.profile.map { meta, f -> f }
        .collectFile(name: 'output_profile_mut_GONO_NGSTAR.csv', keepHeader: true, sort: { it.text })
    rrna23s_collated = RRNA23S.out.profile.map { meta, f -> f }
        .collectFile(name: 'output_profile_GONO_rRNA23S.csv', keepHeader: true, sort: { it.text })
    mlst_collated = MLST_MLST.out.profile.map { meta, f -> f }
        .collectFile(name: 'output_profile_GONO_MLST.csv', keepHeader: true, sort: { it.text })
    ngstar_collated = MLST_NGSTAR.out.profile.map { meta, f -> f }
        .collectFile(name: 'output_profile_GONO_NGSTAR.csv', keepHeader: true, sort: { it.text })
    ngmast_collated = MLST_NGMAST.out.profile.map { meta, f -> f }
        .collectFile(name: 'output_profile_GONO_NGMAST.csv', keepHeader: true, sort: { it.text })

    dna_fasta = MASTERBLASTR_AMR.out.dna.map { meta, f -> f }
        .collectFile(name: 'output_dna.fasta', sort: { it.text }, storeDir: "${params.outdir}/GONO")
    aa_fasta = MASTERBLASTR_AMR.out.aa.map { meta, f -> f }
        .collectFile(name: 'output_aa.fasta', sort: { it.text }, storeDir: "${params.outdir}/GONO")

    // AMR_ALL chain -> gono AMR labware
    LABWARE_GONO_AMR(amr_collated, ngstar_mut_collated, rrna23s_collated, ref_amr)

    // good/bad split for each typing scheme (one MLST_SPLIT per (scheme) item)
    ch_typing = mlst_collated.map   { f -> tuple('GONO', 'MLST',   f) }
        .mix( ngstar_collated.map   { f -> tuple('GONO', 'NGSTAR', f) },
              ngmast_collated.map   { f -> tuple('GONO', 'NGMAST', f) } )
    MLST_SPLIT(ch_typing)

    ch_versions = ch_versions.mix(
        LABWARE_GONO_AMR.out.versions,
        MLST_SPLIT.out.versions.first()
    )

    emit:
    amr      = LABWARE_GONO_AMR.out.labware
    typing   = MLST_SPLIT.out.labware
    versions = ch_versions
}
