process LABWARE_GAS_TOXINS {
    tag "GAS:TOXINS"
    label 'process_single'

    input:
    path profile

    output:
    path "LabWareUpload_GAS_TOXINS.csv", emit: labware
    path "versions.yml",                 emit: versions

    script:
    """
    wade.R labware --org GAS --test TOXINS --profile ${profile} --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r: \$(Rscript -e 'cat(as.character(getRversion()))')
    END_VERSIONS
    """

    stub:
    """
    touch LabWareUpload_GAS_TOXINS.csv
    echo '"${task.process}":' > versions.yml
    """
}
