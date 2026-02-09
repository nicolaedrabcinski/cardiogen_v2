#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { FASTP }                from './modules/fastp'
include { FASTQC }               from './modules/fastqc'
include { BWAMEM_HG38 }          from './modules/bwa_mem'
include { SAMTOOLS_SORT_HG38 }   from './modules/samtools_sort'
include { OCTOPUS_HG38_BWAMEM }  from './modules/octopus'
include { MANTA_HG38_BWAMEM }    from './modules/manta'
include { DELLY_HG38_BWAMEM }    from './modules/delly'

workflow {
    // Устанавливаем значение по умолчанию для min_reads если не задано
    if (!params.min_reads) {
        params.min_reads = 1000
    }
    
    log.info """
    ================================================================
                Cardiogen Illumina WES Data Processing Workflow
    ================================================================
    Input directory : ${params.input_dir}
    Output directory: ${params.outdir}
    File pattern    : ${params.pattern}
    Reference genome: hg38
    Preprocessing  : FastP
    FastQC         : ${params.run_fastqc ? 'Enabled' : 'Disabled'}
    Min reads threshold: ${params.min_reads}
    ================================================================
    """

    fastq_ch = Channel
        .fromFilePairs("${params.input_dir}/${params.pattern}", checkIfExists: true)
        .map { sample_id, files ->
            def meta = [
                id: sample_id, 
                single_end: files.size() == 1,
                qc_tool: 'fastp'
            ]
            tuple(meta, files)
        }

    if (params.run_fastqc) {
        FASTQC(fastq_ch)
    }

    fastp_out = FASTP(fastq_ch)
    
    fastp_reads_ch = fastp_out.map { meta, cleaned_reads, json, html ->
        def read1 = cleaned_reads[0]
        def read2 = cleaned_reads.size() > 1 ? cleaned_reads[1] : null
        tuple(meta, read1, read2)
    }
    
    ref_dir = Channel.value(file(params.hg38_dir))
    bwa_sam_ch = BWAMEM_HG38(fastp_reads_ch, ref_dir)
    bam_bai_ch = SAMTOOLS_SORT_HG38(bwa_sam_ch)
    
    bam_stats_ch = bam_bai_ch
        .map { meta, bam, bai ->
            def read_count = 0
            try {
                def proc = "samtools view -c ${bam}".execute()
                read_count = proc.text.trim().toInteger()
                proc.waitFor()
            } catch (Exception e) {
                log.warn "Cannot count reads for ${meta.id}: ${e.message}"
            }
            
            [meta: meta, bam: bam, bai: bai, read_count: read_count]
        }

    non_empty_bam_bai_ch = bam_stats_ch
        .filter { data -> data.read_count >= params.min_reads }
        .map { data -> 
            log.info "Including sample ${data.meta.id} for variant calling (${data.read_count} reads)"
            tuple(data.meta, data.bam, data.bai) 
        }

    bam_stats_ch
        .filter { data -> data.read_count < params.min_reads && data.read_count > 0 }
        .subscribe { data ->
            log.warn "Excluding sample ${data.meta.id} from variant calling: low coverage (${data.read_count} reads < ${params.min_reads})"
        }

    bam_stats_ch
        .filter { data -> data.read_count == 0 }
        .subscribe { data ->
            log.warn "Excluding sample ${data.meta.id} from variant calling: empty BAM (0 reads)"
        }

    ref_fasta = Channel.value(file("${params.hg38_dir}/hg38.analysisSet.fa"))
    ref_fai = Channel.value(file("${params.hg38_dir}/hg38.analysisSet.fa.fai"))

    OCTOPUS_HG38_BWAMEM(non_empty_bam_bai_ch, ref_fasta, ref_fai)
    MANTA_HG38_BWAMEM(non_empty_bam_bai_ch, ref_fasta, ref_fai)
    DELLY_HG38_BWAMEM(non_empty_bam_bai_ch, ref_fasta, ref_fai)
}

workflow.onComplete {
    log.info """
    ================================================================
                    Pipeline completed successfully!
    ================================================================
    Results saved in: ${params.outdir ?: 'results'}/
    ================================================================
    """
}