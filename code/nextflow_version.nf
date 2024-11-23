#!/usr/bin/env nextflow

process FASTQC {
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("fastqc_result")

    script:
    """
    mkdir -p fastqc_result
    fastqc -o fastqc_result ${reads}
    """
}

process MULTIQC {
    publishDir params.output_dir, mode:'copy'

    input:
    tuple val(sample_id), path(fastqc)

    output:
    path "multiQC_report.html"

    script:
    """
    multiqc --filename multiQC_report $fastqc
    """
}

process STARsolo {
    input:
    tuple val(sample_id), path(reads)
    path genome_dir
    path whitelist
    path output_dir
    val threads

    output:
    tuple val(sample_id), path("tmp/STARsolo")

    script:
    """
    mkdir -p "tmp/STARsolo/"
    STAR --runThreadN ${threads} \
         --genomeDir ${genome_dir} \
         --readFilesIn ${reads} \
         --soloType CB_UMI_Simple \
         --soloCBwhitelist ${whitelist} \
         --soloUMIlen 12 \
         --readFilesCommand zcat \
         --outFileNamePrefix tmp/STARsolo/
    """
}

process CompressFiles {
    input:
    tuple val(sample_id), path(raw_file)

    output:
    tuple val(sample_id), path("${raw_file}/Solo.out/Gene/raw/")

    script:
    """
    gzip ${raw_file}/Solo.out/Gene/raw/*
    """
}

process seurat_obj {
    input:
    tuple val(sample_id), path(data_dir)
    path script
    path output_dir

    output:
    tuple val(sample_id), path("seurat_obj.rds")

    script:
    """
    mkdir -p ${output_dir}
    Rscript $script ${data_dir} "."
    """
}

process Draw_UMAP {
    publishDir params.output_dir, mode: 'copy'

    input:
    tuple val(sample_id), path(seurat_obj)
    path script
    path output_dir

    output:
    tuple val(sample_id), path("umap_plot.png")

    script:
    """
    mkdir -p ${output_dir}
    Rscript $script ${seurat_obj} "."
    """
}

process Draw_PCA {
    publishDir params.output_dir, mode: 'copy'

    input:
    tuple val(sample_id), path(seurat_obj)
    path script
    path output_dir

    output:
    tuple val(sample_id), path("pca_plot.png")

    script:
    """
    mkdir -p ${output_dir}
    Rscript $script ${seurat_obj} "."
    """
}

workflow {
    fastq_files_ch = Channel
        .fromFilePairs("${params.fastq_dir}/*_{1,2}.fastq.gz", checkIfExists: true)
        .map { pair -> [pair[0], pair[1].reverse()] }

    FASTQC_ch = FASTQC(fastq_files_ch)

    MULTIQC(FASTQC_ch)

    STARsolo_ch = STARsolo(fastq_files_ch, params.genome_dir, params.whitelist, params.output_dir, params.threads)
    // STARsolo_ch.view()

    CompressFiles_ch = CompressFiles(STARsolo_ch)

    seurat_obj_ch = seurat_obj(CompressFiles_ch, params.seurat_obj_script, params.output_dir)

    Draw_UMAP(seurat_obj_ch, params.UMAP_script, params.output_dir)

    Draw_PCA(seurat_obj_ch, params.PCA_script, params.output_dir)
}

