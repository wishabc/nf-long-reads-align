#!/usr/bin/env nextflow


def set_key_for_group_tuple(ch) {
  ch.groupTuple()
  .map{ it -> tuple(groupKey(it[0], it[1].size()), *it[1..(it.size()-1)]) }
  .transpose()
}


process align_reads {
    publishDir "${params.outdir}/minimap_logs", pattern: "log.*"
    tag "${cell_line}:${name}"
    conda params.conda

    input:
        tuple val(cell_line), path(fastq)
    output:
        tuple val(cell_line), path(name)
    
    script:
    name = "${fastq.simpleName}.sam"
    log_files = "log.${fastq.simpleName}.txt"
    """
    minimap2 ${params.genome} ${fastq} -a -z 600,200 -x map-ont --MD -Y -o ${name} &>${log_files}
    """
}

process merge_files {
    scratch true
    publishDir "${params.outdir}/alignments"
    tag "${cell_line}"
    conda params.conda

    input:
        tuple val(cell_line), path(bams)
    output:
        tuple val(cell_line), path(name), path("${name}.bai")
    
    script:
    name = "${cell_line}.bam"
    """
    samtools merge -u merged.bam ${bams}
    samtools sort -O bam -o ${name} merged.bam
    samtools index ${name}
    """
}

workflow {
    fasta_ch = set_key_for_group_tuple(
        Channel.fromPath(params.samples_file)
        .splitCsv(header:true, sep:'\t')
        .map(
            row -> tuple(row.cell_line, file(row.fastq_file))
            )
        )
    
    al_reads = align_reads(fasta_ch)
    merge_files(al_reads.groupTuple())
}