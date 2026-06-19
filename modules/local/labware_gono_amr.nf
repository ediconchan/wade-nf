process LABWARE_GONO_AMR {
    tag "GONO:AMR"
    label 'process_single'

    input:
    path amr_profile
    path ngstar_profile
    path rrna23s_profile
    path refdir

    output:
    path "LabWareUpload_GONO_AMR.csv",      emit: labware
    path "LabWareUpload_GONO_AMR_good.csv", emit: good
    path "LabWareUpload_GONO_AMR_bad.csv",  emit: bad
    path "Elevated_MICs.csv",               emit: elevated
    path "versions.yml",                    emit: versions

    script:
    """
    wade.R labware --org GONO --test AMR \\
        --profile ${amr_profile} \\
        --profile2 ${ngstar_profile} \\
        --profile3 ${rrna23s_profile} \\
        --refdir ${refdir} \\
        --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r: \$(Rscript -e 'cat(as.character(getRversion()))')
    END_VERSIONS
    """

    stub:
    """
    touch LabWareUpload_GONO_AMR.csv LabWareUpload_GONO_AMR_good.csv LabWareUpload_GONO_AMR_bad.csv Elevated_MICs.csv
    echo '"${task.process}":' > versions.yml
    """
}
