process DELLY_HG38_BWAMEM {
    tag "$meta.id - $meta.qc_tool - hg38 - bwamem"
    label 'process_high'

    publishDir "${params.outdir}/06_variant_calling/delly/bwamem/${meta.qc_tool}/hg38", mode: 'copy'

    container 'cardiogen-wes-tools:latest'

    input:
    tuple val(meta), path(bam), path(bai)
    path reference
    path reference_index

    output:
    tuple val(meta), path("*.bcf"),     emit: vcf
    tuple val(meta), path("*.bcf.csi"), emit: tbi
    path "versions.yml",                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.id}_${meta.qc_tool}_hg38_bwamem"
    
    """
    # Run DELLY variant calling
    delly call \\
        -g ${reference} \\
        -o ${prefix}.bcf \\
        $args \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        delly: \$(delly 2>&1 | grep "Version:" | sed 's/^.*Version: //; s/).*\$//')
    END_VERSIONS
    """
}