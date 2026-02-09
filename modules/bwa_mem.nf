process BWAMEM_HG38 {
    tag "$meta.id"
    
    input:
    tuple val(meta), path(read1), path(read2)
    path ref_dir
    
    output:
    tuple val(meta), path("${meta.id}.sam")
    
    script:
    def read_args = read2 ? "${read1} ${read2}" : "${read1}"
    def ref_path = "${ref_dir}/hg38.analysisSet.fa"
    
    """
    bwa mem -t $task.cpus \\
        -R "@RG\\\\tID:${meta.id}\\\\tSM:${meta.id}\\\\tPL:ILLUMINA\\\\tLB:${meta.id}" \\
        "$ref_path" \\
        $read_args \\
        > ${meta.id}.sam
    """
}