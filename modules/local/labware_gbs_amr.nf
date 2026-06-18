process LABWARE_GBS_AMR {
    tag "GBS:AMR"
    label 'process_single'

    input:
    path profile

    output:
    path "LabWareUpload_GBS_AMR.csv",     emit: labware
    path "LabWareUpload_GBS_AMR_bad.csv", emit: bad
    path "versions.yml",                  emit: versions

    script:
    """
    wade.R labware --org GBS --test AMR --profile ${profile} --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r: \$(Rscript -e 'cat(as.character(getRversion()))')
    END_VERSIONS
    """

    stub:
    """
    touch LabWareUpload_GBS_AMR.csv LabWareUpload_GBS_AMR_bad.csv
    echo '"${task.process}":' > versions.yml
    """
}
