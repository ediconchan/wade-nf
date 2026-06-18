process LABWARE_PNEUMO_AMR {
    tag "PNEUMO:AMR"
    label 'process_single'

    input:
    path amr_profile
    path rrna23s_profile
    path refdir

    output:
    path "LabWareUpload_PNEUMO_AMR.csv",      emit: labware
    path "LabWareUpload_PNEUMO_AMR_good.csv", emit: good
    path "LabWareUpload_PNEUMO_AMR_bad.csv",  emit: bad
    path "versions.yml",                      emit: versions

    script:
    """
    wade.R labware --org PNEUMO --test AMR \\
        --profile ${amr_profile} \\
        --profile2 ${rrna23s_profile} \\
        --refdir ${refdir} \\
        --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r: \$(Rscript -e 'cat(as.character(getRversion()))')
    END_VERSIONS
    """

    stub:
    """
    touch LabWareUpload_PNEUMO_AMR.csv LabWareUpload_PNEUMO_AMR_good.csv LabWareUpload_PNEUMO_AMR_bad.csv
    echo '"${task.process}":' > versions.yml
    """
}
