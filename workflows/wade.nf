//
// WADE main workflow: parse samplesheet, branch by organism, run per-org routines
//
include { GAS } from '../subworkflows/local/gas.nf'
include { GBS } from '../subworkflows/local/gbs.nf'

workflow WADE {

    take:
    ch_input   // path to samplesheet csv
    wade_data  // path to wade-data base dir

    main:
    ch_versions = Channel.empty()

    //
    // Parse samplesheet: sample,organism,assembly,vcf,variable
    // `organism` falls back to params.organism when blank. Relative assembly/vcf
    // paths are resolved against the samplesheet's own directory so bundled test
    // data works regardless of launch directory.
    //
    def sheet_dir = file(ch_input).parent
    def resolve_path = { p ->
        // absolute paths / URIs pass through; relative paths resolve against the
        // samplesheet's directory (file() would otherwise resolve to launchDir)
        (p ==~ /^(\/|[a-z0-9]+:\/\/).*/) ? file(p) : file("${sheet_dir}/${p}")
    }

    ch_rows = Channel.fromPath(ch_input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def org = (row.organism && row.organism.trim()) ? row.organism.trim().toUpperCase() : (params.organism ? params.organism.toUpperCase() : null)
            if (!org) {
                error "Sample '${row.sample}': no organism in samplesheet and --organism not set"
            }
            def meta = [ id: row.sample, org: org, variable: (row.variable ?: '') ]
            def assembly = resolve_path(row.assembly)
            if (!assembly.exists()) error "Sample '${row.sample}': assembly not found: ${row.assembly}"
            def vcf = (row.vcf && row.vcf.trim()) ? resolve_path(row.vcf.trim()) : []
            [ meta, assembly, vcf ]
        }

    // contig channel [meta, assembly]
    ch_contigs = ch_rows.map { meta, assembly, vcf -> [ meta, assembly ] }

    //
    // Branch by organism
    //
    ch_branched = ch_contigs.branch { meta, assembly ->
        gas:    meta.org == 'GAS'
        gbs:    meta.org == 'GBS'
        gono:   meta.org == 'GONO'
        pneumo: meta.org == 'PNEUMO'
        other:  true
    }

    ch_branched.other.map { meta, assembly ->
        log.warn "Organism '${meta.org}' (sample ${meta.id}) is not supported; skipping"
    }

    //
    // Per-organism routines
    //
    GAS(ch_branched.gas, wade_data)
    GBS(ch_branched.gbs, wade_data)
    ch_versions = ch_versions.mix(GAS.out.versions, GBS.out.versions)

    //
    // GONO / PNEUMO subworkflows are wired the same way once their LabWare
    // formatters are ported (SPEC.md migration plan phase 3-4).
    //

    //
    // Collect software versions
    //
    ch_versions
        .unique()
        .collectFile(name: 'software_versions.yml', storeDir: "${params.outdir}/pipeline_info")

    emit:
    versions = ch_versions
}
