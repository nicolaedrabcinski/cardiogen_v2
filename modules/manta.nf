process MANTA_HG38_BWAMEM {
    tag "$meta.id - $meta.qc_tool - hg38 - bwamem"
    label 'process_high'

    publishDir "${params.outdir}/06_variant_calling/manta/bwamem/${meta.qc_tool}/hg38", mode: 'copy'

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
    # Configure Manta
    configManta.py \\
        --bam ${bam} \\
        --referenceFasta ${reference} \\
        --runDir manta_work \\
        $args

    # Run Manta
    manta_work/runWorkflow.py -m local -j ${task.cpus}

    # Copy results
    cp manta_work/results/variants/diploidSV.vcf.gz ${prefix}.vcf.gz
    cp manta_work/results/variants/diploidSV.vcf.gz.tbi ${prefix}.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        manta: \$(configManta.py --version 2>&1 | grep -oP 'Manta workflow version: \\K[0-9.]+')
    END_VERSIONS
    """
}