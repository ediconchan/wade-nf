process EMM {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(contig)
    path refdir

    output:
    tuple val(meta), path("output_profile_emm.csv"), emit: profile
    path "versions.yml",                              emit: versions

    script:
    def prefix   = task.ext.prefix ?: meta.id
    def variable = meta.variable ?: ''
    """
    wade.R emm \\
        --org ${meta.org} \\
        --sample ${prefix} \\
        --variable '${variable}' \\
        --contig ${contig} \\
        --refdir ${refdir} \\
        --outdir . \\
        --tempdir wade_tmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | sed -n 's/^blastn: //p')
        r: \$(Rscript -e 'cat(as.character(getRversion()))')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    echo "Sample,Type,Subtype_Rep,Subtype,bp_id,Comments" > output_profile_emm.csv
    echo "${prefix},emm1,emm1.0,emm1.0,180," >> output_profile_emm.csv
    echo '"${task.process}":' > versions.yml
    """
}
