#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    wade-nf : Nextflow port of phac-nml/wade (molecular typing of GAS/GBS/PNEUMO/GONO)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

include { WADE } from './workflows/wade.nf'

workflow {
    if (!params.input) {
        error "You must provide a samplesheet with --input <samplesheet.csv>"
    }
    if (!params.wade_data) {
        error "You must provide the reference data dir with --wade_data <dir> (the phac-nml/wade-data layout)"
    }

    WADE(params.input, params.wade_data)
}
