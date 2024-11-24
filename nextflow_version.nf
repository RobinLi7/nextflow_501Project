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

process download {
    input:
    path ref_dir

    output:
    path("ref")

    script:
    """
    mkdir -p "ref/"
    wget -O ${ref_dir}/CellRanger/3M-february-2018.txt.gz https://github.com/noamteyssier/10x_whitelist_mirror/raw/main/3M-february-2018.txt.gz
    gunzip -f ${ref_dir}/CellRanger/3M-february-2018.txt.gz
    wget -O ${ref_dir}/chr1.fa.gz http://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr1.fa.gz
    wget -O ${ref_dir}/chr2.fa.gz http://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr2.fa.gz
    wget -O ${ref_dir}/chr3.fa.gz http://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr3.fa.gz
    wget -O ${ref_dir}/gencode.v38.annotation.gtf.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz 
    gunzip -f ${ref_dir}/chr*.fa.gz
    gunzip -f ${ref_dir}/*.gtf.gz
    grep -P '^chr1\t|^chr2\t|^chr3\t' ${ref_dir}/gencode.v38.annotation.gtf > ${ref_dir}/chr1_chr2_chr3.gtf
    """
}

process ref_prep {
    // publishDir params.genome_dir, mode:'copy'

    input:
    path ref_path
    path genome_dir
    path gtf_file
    val threads

    output:
    path("ref/STAR_genomeDir")

    script:
    """
    mkdir -p "ref/STAR_genomeDir/"
    STAR --runThreadN ${threads} \
     --runMode genomeGenerate \
     --genomeDir ref/STAR_genomeDir/ \
     --genomeFastaFiles ${ref_path}/chr1.fa ${ref_path}/chr2.fa ${ref_path}/chr3.fa \
     --sjdbGTFfile ${gtf_file} \
     --sjdbOverhang 90 \
     --genomeSAindexNbases 13/
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
         --soloCBwhitelist ${whitelist} \
         --soloUMIlen 12 \
         --readFilesCommand zcat \
         --soloType CB_UMI_Simple \
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

    ref_raw_path = download(params.ref_dir)

    ref_path = ref_prep(ref_raw_path, params.genome_dir, params.GTF_FILE, params.threads)

    STARsolo_ch = STARsolo(fastq_files_ch, ref_path, params.whitelist, params.output_dir, params.threads)
    // STARsolo_ch = STARsolo(fastq_files_ch, "/projects/rli_prj/large_files/STAR_genomeDir", params.whitelist, params.output_dir, params.threads)

    CompressFiles_ch = CompressFiles(STARsolo_ch)

    seurat_obj_ch = seurat_obj(CompressFiles_ch, params.seurat_obj_script, params.output_dir)

    Draw_UMAP(seurat_obj_ch, params.UMAP_script, params.output_dir)

    Draw_PCA(seurat_obj_ch, params.PCA_script, params.output_dir)
}

