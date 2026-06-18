process MLST_SPLIT {
    tag "${org}:${test}"
    label 'process_single'

    input:
    tuple val(org), val(test), path(profile)

    output:
    path "LabWareUpload_${org}_${test}.csv",      emit: labware
    path "LabWareUpload_${org}_${test}_good.csv", emit: good
    path "LabWareUpload_${org}_${test}_bad.csv",  emit: bad
    path "versions.yml",                          emit: versions

    script:
    """
    wade.R mlst_split --org ${org} --test ${test} --profile ${profile} --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r: \$(Rscript -e 'cat(as.character(getRversion()))')
    END_VERSIONS
    """

    stub:
    """
    touch LabWareUpload_${org}_${test}.csv LabWareUpload_${org}_${test}_good.csv LabWareUpload_${org}_${test}_bad.csv
    echo '"${task.process}":' > versions.yml
    """
}
