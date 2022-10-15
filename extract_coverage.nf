#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { set_key_for_group_tuple } from "./aligning"


process split_region {
    tag "${region}"
    scratch true
    conda params.conda

    input:
        val region
    output:
        path "${bedfile}*"

    script:
    bedfile = region + '.atomic.bed.'
    region_bed = region.replaceAll('_', '\t')
    """
    echo "${region_bed}" > temp.bed
    bedops --chop "${params.slice_size}" temp.bed | bedops -n 1 - "${params.blacklist}" | head -n -1 | split -l 500 - ${bedfile}
    """
}

process mpileup {
    tag "${id}:${region}"
    conda params.conda
    scratch true

    input:
        tuple val(id), val(bam_file), val(bam_file_index), path(bedfile)
    output:
        tuple val(id), path(count_list)

    script:
    region = "${bedfile.simpleName}"
    count_list = "${id}.${region}.counts.bed"
    """
    python3 ${projectDir}/bin/count_tags.py "${bam_file}" < "${bedfile}" > tags_counts.txt
    paste "${bedfile}" tags_counts.txt > ${count_list} 
    """
}

process sort {
    tag "${id}"
    publishDir params.outdir
    scratch true
    conda params.conda

    input:
        tuple val(id), path(counts)
    output:
        path result_bed

    script:
    result_bed = "${id}.counts_by_chunks.bed"
    counts_list = counts.split(' ').join('\n')
    """
    echo "${counts_list}" > filenames.txt
    while read line; do
        cat \$line >> collected_file.bed
    done < filenames.txt

    sort -k 1,1 -k2,2n collected_file.bed > ${result_bed}
    """
}

workflow {
    bam_files = Channel.fromPath(params.bams_list)
        .splitText()
        .map(it -> it.strip())
        .map(
            it -> tuple(file(it).simpleName, file(it), file("${it}.*ai"))
        )
    
    input = Channel.fromPath(params.nanosv_regions).splitText()
        .map(it -> it.trim().replaceAll('\t', '_'))

    regions = split_region(input) | flatten
    pileup = mpileup(set_key_for_group_tuple(bam_files.combine(regions)))
    
    sort(pileup.groupTuple())
}