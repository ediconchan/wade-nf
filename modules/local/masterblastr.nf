process MASTERBLASTR {
    tag "${meta.id}:${test}"
    label 'process_single'

    input:
    tuple val(meta), path(contig)
    path refdir
    val  test

    output:
    tuple val(meta), path("output_profile_*.csv"), emit: profile
    tuple val(meta), path("output_dna.fasta"),                  optional: true, emit: dna
    tuple val(meta), path("output_aa.fasta"),                   optional: true, emit: aa
    tuple val(meta), path("output_dna_notfound*.fasta"),        optional: true, emit: notfound
    path "versions.yml",                            emit: versions

    script:
    def prefix   = task.ext.prefix ?: meta.id
    def variable = meta.variable ?: ''
    """
    wade.R masterblastr \\
        --org ${meta.org} \\
        --test ${test} \\
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
    echo "SampleNo,Variable,stub_result,SampleProfile" > output_profile_${meta.org}_${test}.csv
    echo "${prefix},,POS,Susceptible" >> output_profile_${meta.org}_${test}.csv
    touch output_dna.fasta output_aa.fasta
    echo '"${task.process}":' > versions.yml
    """
}
