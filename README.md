# nextflow_501Project

## Background and Rationale

This pipeline is designed for single-cell RNA sequencing (scRNA-seq) data analysis. It uses **FastQC**, **MultiQC**, **STARsolo**, and **Seurat** to process, analyze, and visualize scRNA-seq data. Key features include quality control, genome generation, mapping, and downstream clustering analysis with visualization.

### Aims:
- Perform quality control of raw sequencing reads.
- Generate genome references and map reads using STARsolo.
- Create and analyze Seurat objects for clustering and visualization.

### Dependencies:
- git
- Nextflow (`=21.10.0`)
- Singularity (`=3.5.2-1.1.el7`for containerized execution)
- Software dependencies included in containers:
  - **biocontainers/fastqc**
  - **ewels/multiqc**
  - **quay.io/biocontainers/star**
  - **satijalab/seurat**

### DAG
![DAG](DAG.png)
---

## Usage

### Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/RobinLi7/nextflow_501Project.git
   ```
   
2. Make sure Nextflow (`=21.10.0`) is installed.\
   [Installation guide](https://www.nextflow.io/docs/latest/install.html)

4. Make sure Singularity (`=3.5.2-1.1.el7`) is installed and configured. If you are using a BCGSC computer, there is no need for installation or configuration.\
   [Installation guide](https://docs.sylabs.io/guides/3.0/user-guide/installation.html)
   
6. Download the reference genome files:\
   There is a `download` process in my Nextflow pipeline that will handle the downloads automatically.

---

### Running the Pipeline
Go to the project directory
```bash
cd nextflow_501Project/
```
And run
```bash
nextflow run nextflow_version.nf 
```

### Parameters

| Parameter            | Description                                               | Default Value                              |
|-----------------------|-----------------------------------------------------------|-------------------------------------------|
| `fastq_dir`          | Directory containing the input FASTQ files.                | `${projectDir}/data`                      |
| `output_dir`         | Directory where the output files will be stored.           | `${projectDir}/result`                    |
| `ref_dir`            | Directory for storing reference files.                     | `${projectDir}/ref`                       |
| `GTF_FILE`           | Path to the GTF file for gene annotations.                 | `${projectDir}/ref/chr1_chr2_chr3.gtf`    |
| `genome_dir`         | Directory containing the STAR genome index.                | `${projectDir}/ref/STAR_genomeDir`        |
| `threads`            | Number of threads for parallel execution.                  | `8`                                       |
| `whitelist`          | Path to the whitelist file used for cell barcodes.         | `${projectDir}/ref/CellRanger/3M-february-2018.txt` |
| `seurat_obj_script`  | Path to the R script for creating Seurat objects.          | `${projectDir}/code/seurat_obj.R`         |
| `UMAP_script`        | Path to the R script for generating UMAP plots.            | `${projectDir}/code/draw_UMAP.R`          |
| `PCA_script`         | Path to the R script for generating PCA plots.             | `${projectDir}/code/draw_PCA.R`           |

---

### Updated Input Explanation

The pipeline processes **10x Genomics scRNA-seq data** and requires the following input data:

### 1. **FASTQ Files**
- The input FASTQ files are generated from 10x Genomics platforms and follow this standard naming convention:
  ```
  pbmc_1k_v3_S1_L001_R1_001.fastq.gz
  pbmc_1k_v3_S1_L001_R2_001.fastq.gz
  ```
  - **R1**: Contains the cell barcodes and unique molecular identifiers (UMIs).
  - **R2**: Contains the actual RNA sequences for downstream analysis.

- For demonstration purposes, this pipeline uses a **subsample** of the original data. The subsampling is performed to reduce computational resource usage while testing the pipeline. 

- **Subsampling Command**: The following `seqtk` command was used to randomly sample 2,000,000 reads from the original FASTQ file:
  ```bash
  seqtk sample -s777 /projects/rli_prj/CPSC501/project/data/pbmc_1k_v3_S1_L001_R1_001.fastq.gz 2000000 > /projects/rli_prj/CPSC501/project/data/genome_1.fastq
  ```
  - **`-s777`**: Sets the random seed for reproducibility.
  - **`2000000`**: Specifies the number of reads to sample.
  - The resulting file `genome_1.fastq` is included in the `data` folder for this pipeline.

### 2. **Genome Files**
The pipeline uses a subset of the human genome as a reference, specifically the chromosomes `chr1`, `chr2`, and `chr3`. These files include:
- **FASTA Files**: Contain the nucleotide sequences for the reference genome. The required files are:
  - `chr1.fa`: Sequence for chromosome 1.
  - `chr2.fa`: Sequence for chromosome 2.
  - `chr3.fa`: Sequence for chromosome 3.
- **Annotation File (GTF)**: Provides functional annotations of genes, including their locations on the genome and exon structures.
  - File used: `gencode.v38.annotation.gtf`.
  - This file is filtered to include annotations only for `chr1`, `chr2`, and `chr3`.

### 3. **Barcode Whitelist**
- A whitelist file is required for cell barcode validation during alignment. 
- File: `3M-february-2018.txt`, provided by the Cell Ranger dataset.
- It ensures that only valid cell barcodes are used in the downstream analysis.

### 4. **Optional Custom Inputs**
- If users wish to expand the analysis to additional chromosomes or use a different reference genome, they need to replace the provided genome and annotation files with their own.
- Modify the `GTF_FILE` and `genome_dir` parameters in `nextflow.config` to reflect the new files.

### Important Notes:
- All genome and annotation files are downloaded and preprocessed automatically by the pipeline during execution.
- Ensure that the input files are correctly named and formatted to avoid errors during processing.

---

## Output

- **Quality Control**:
  - Per-sample FastQC reports (`fastqc_result`).
  - MultiQC summary (`multiQC_report.html`).

- **Genome Files**:
  - STAR genome index files in `ref/STAR_genomeDir`.

- **Mapping Results**:
  - STARsolo outputs in `tmp/STARsolo`.

- **Seurat Objects**:
  - Processed object: `seurat_obj.rds`.

- **Visualizations**:
  - UMAP plot (`umap_plot.png`).
  - PCA plot (`pca_plot.png`).

---

## Notes

- Ensure proper Singularity setup for containerized execution.
- For any issues, please contact contributors listed in `nextflow.config`.
