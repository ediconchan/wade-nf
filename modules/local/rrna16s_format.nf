process RRNA16S_FORMAT {
    tag "${org}"
    label 'process_single'

    input:
    tuple val(org), path(profile)

    output:
    path "LabWareUpload_*_rRNA16S.csv", emit: labware
    path "versions.yml",                emit: versions

    script:
    """
    wade.R rrna16s --org ${org} --profile ${profile} --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r: \$(Rscript -e 'cat(as.character(getRversion()))')
    END_VERSIONS
    """

    stub:
    """
    echo "SampleNo,Allele,Identification,Percent_Match" > LabWareUpload_${org}_rRNA16S.csv
    echo '"${task.process}":' > versions.yml
    """
}
