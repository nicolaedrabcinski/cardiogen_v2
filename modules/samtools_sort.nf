process SAMTOOLS_SORT_HG38 {
    tag "$meta.id"
    
    input:
    tuple val(meta), path(sam_file)
    
    output:
    tuple val(meta), path("${meta.id}_hg38.bam"), path("${meta.id}_hg38.bam.bai")
    
    publishDir "${params.outdir}/bam", mode: 'copy'  // копировать в results/bam

    script:
    """
    samtools sort -@ $task.cpus -o ${meta.id}_hg38.bam ${sam_file}
    samtools index -@ $task.cpus ${meta.id}_hg38.bam
    """
}