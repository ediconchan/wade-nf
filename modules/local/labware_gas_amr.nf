process LABWARE_GAS_AMR {
    tag "GAS:AMR"
    label 'process_single'

    input:
    path profile

    output:
    path "LabWareUpload_GAS_AMR.csv",     emit: labware
    path "LabWareUpload_GAS_AMR_bad.csv", emit: bad
    path "versions.yml",                  emit: versions

    script:
    """
    wade.R labware --org GAS --test AMR --profile ${profile} --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r: \$(Rscript -e 'cat(as.character(getRversion()))')
    END_VERSIONS
    """

    stub:
    """
    touch LabWareUpload_GAS_AMR.csv LabWareUpload_GAS_AMR_bad.csv
    echo '"${task.process}":' > versions.yml
    """
}
