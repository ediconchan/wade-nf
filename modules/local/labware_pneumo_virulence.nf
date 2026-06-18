process LABWARE_PNEUMO_VIRULENCE {
    tag "PNEUMO:VIRULENCE"
    label 'process_single'

    input:
    path profile

    output:
    path "LabWareUpload_PNEUMO_VIRULENCE.csv", emit: labware
    path "versions.yml",                       emit: versions

    script:
    """
    wade.R labware --org PNEUMO --test VIRULENCE --profile ${profile} --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r: \$(Rscript -e 'cat(as.character(getRversion()))')
    END_VERSIONS
    """

    stub:
    """
    touch LabWareUpload_PNEUMO_VIRULENCE.csv
    echo '"${task.process}":' > versions.yml
    """
}
