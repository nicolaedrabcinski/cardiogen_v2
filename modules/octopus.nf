process OCTOPUS_HG38_BWAMEM {
    tag "$meta.id - $meta.qc_tool - hg38 - bwamem"
    label 'process_high'
    
    publishDir "${params.outdir}/06_variant_calling/octopus/bwamem/${meta.qc_tool}/hg38", mode: 'copy'

    container 'cardiogen-wes-tools:latest'

    input:
    tuple val(meta), path(bam), path(bai)
    path reference
    path reference_index

    output:
    tuple val(meta), path("*.vcf.gz"),     emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
    path "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.id}_${meta.qc_tool}_hg38_bwamem"
    
    """
    # Run Octopus variant calling (compatible with v0.7.4)
    octopus \\
        -R ${reference} \\
        -I ${bam} \\
        -o ${prefix}.vcf.gz \\
        --threads ${task.cpus} \\
        --max-haplotypes 200 \\
        --downsample-above 1000 \\
        --downsample-target 500 \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        octopus: \$(octopus --version 2>&1 | sed 's/octopus //g')
    END_VERSIONS
    """
}