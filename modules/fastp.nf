process FASTP {
    tag "$meta.id"
    
    input:
    tuple val(meta), path(reads)
    
    output:
    tuple val(meta), path("*.fastq.gz"), path("*.json"), path("*.html")
    
    script:
    if (meta.single_end) {
        """
        fastp -i ${reads} -o ${meta.id}_fastp_1.fastq.gz \\
            --json ${meta.id}_fastp.json \\
            --html ${meta.id}_fastp.html \\
            $params.fastp_args
        """
    } else {
        """
        fastp -i ${reads[0]} -I ${reads[1]} \\
            -o ${meta.id}_fastp_1.fastq.gz -O ${meta.id}_fastp_2.fastq.gz \\
            --json ${meta.id}_fastp.json \\
            --html ${meta.id}_fastp.html \\
            $params.fastp_args
        """
    }
}