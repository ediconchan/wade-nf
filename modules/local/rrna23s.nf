process RRNA23S {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("output_profile_*_rRNA23S.csv"), emit: profile
    path "versions.yml",                                   emit: versions

    script:
    def prefix = task.ext.prefix ?: meta.id
    def vcf_arg = vcf ? "--vcf ${vcf}" : ""
    """
    wade.R rrna23s \\
        --org ${meta.org} \\
        --sample ${prefix} \\
        ${vcf_arg} \\
        --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r: \$(Rscript -e 'cat(as.character(getRversion()))')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    echo "SampleNo,A2059G,C2611T" > output_profile_${meta.org}_rRNA23S.csv
    echo "${prefix},0,0" >> output_profile_${meta.org}_rRNA23S.csv
    echo '"${task.process}":' > versions.yml
    """
}
