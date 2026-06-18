process SEROTYPE {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(contig)
    path refdir

    output:
    tuple val(meta), path("LabWareUpload_*_SEROTYPE.csv"), emit: profile
    path "versions.yml",                                   emit: versions

    script:
    def prefix = task.ext.prefix ?: meta.id
    """
    wade.R serotype \\
        --org ${meta.org} \\
        --sample ${prefix} \\
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
    echo "SampleNo,Serotype,Profile" > LabWareUpload_${meta.org}_SEROTYPE.csv
    echo "${prefix},III,stub" >> LabWareUpload_${meta.org}_SEROTYPE.csv
    echo '"${task.process}":' > versions.yml
    """
}
