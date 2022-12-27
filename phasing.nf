#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.wh_conda = "/home/sabramov/miniconda3/envs/whatshap"

process extract_genotypes {
    conda params.wh_conda
    tag "${indiv_id}"

    input:
        val indiv_id
    output:
        tuple val(indiv_id), path(name), path("${name}.csi")

    script:
    name = "${indiv_id}.vcf.gz"
    """
    bcftools view -s ${indiv_id} \
        ${params.genotypes_file} | bcftools view -i 'GT == "0/1"' -O z > ${name}
    bcftools index ${name}
    """
}

process whatshap_phase {
    conda params.wh_conda
    tag "${indiv_id}"
    publishDir "${params.outdir}/phased_vcfs", pattern: "${name}*"
    scratch true

    input:
        tuple val(indiv_id), path(bam_file), path(bam_file_index), path(genotype), path(genotype_index)

    output:
        tuple val(indiv_id), path(bam_file), path(bam_file_index), path(name), path("${name}.csi")

    script:
    name = "${indiv_id}.phased.vcf.gz"
    """
    whatshap phase -o phased.vcf --ignore-read-groups \
        --reference=${params.genome_fasta} \
        ${genotype} ${bam_file} 
    bcftools view -O z phased.vcf > ${name}
    bcftools index ${name}
    """
}

process whatshap_haplotag {
    conda params.wh_conda
    tag "${indiv_id}"
    publishDir: "${params.outdir}/phased_bam_files", pattern: "${name}*"
    publishDir: "${params.outdir}/haplotag_stats", pattern: "${stats}"

    input:
        tuple val(indiv_id), path(bam_file), path(bam_file_index), path(phased_vcf), path(phased_vcf_index)

    output:
        tuple val(indiv_id), path(name), path("${name}.bai"), emit: bam
        tuple val(indiv_id), path(stats), emit: stats

    script:
    stats = "${indiv_id}.haplotypes.tsv"
    name = "${indiv_id}.haplotaged.bam"
    """
    whatshap haplotag --output-haplotag-list ${stats} \
        -o ${name} --reference ${params.genome_fasta} \
        --ignore-read-groups --skip-missing-contigs \
        --tag-supplementary \
        ${phased_vcf} ${bam_file}
    """
}


workflow {
    bam_files = Channel.fromPath(params.samples_file)
        .splitCsv(header:true, sep:'\t')
        .map(row -> tuple(row.indiv_id, file(row.bam_file), file("${row.bam_file}.*ai")))

    genotypes = extract_genotypes(bam_files.map(it -> it[0]))
    haplotag = whatshap_phase(bam_files.join(genotypes)) | whatshap_haplotag
}