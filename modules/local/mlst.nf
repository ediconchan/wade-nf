process MLST {
    tag "${meta.id}:${scheme}"
    label 'process_single'

    input:
    tuple val(meta), path(contig)
    path refdir
    val  scheme

    output:
    tuple val(meta), path("output_profile_*.csv"), emit: profile
    path "versions.yml",                            emit: versions

    script:
    def prefix   = task.ext.prefix ?: meta.id
    def variable = meta.variable ?: ''
    """
    wade.R mlst \\
        --org ${meta.org} \\
        --scheme ${scheme} \\
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
    echo "SampleNo,ST" > output_profile_${meta.org}_${scheme}.csv
    echo "${prefix},1" >> output_profile_${meta.org}_${scheme}.csv
    echo '"${task.process}":' > versions.yml
    """
}
