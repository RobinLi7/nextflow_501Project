# nextflow_501Project

## Background and Rationale

This pipeline is designed for single-cell RNA sequencing (scRNA-seq) data analysis. It uses **FastQC**, **MultiQC**, **STARsolo**, and **Seurat** to process, analyze, and visualize scRNA-seq data. Key features include quality control, genome generation, mapping, and downstream clustering analysis with visualization.

### Aims:
- Perform quality control of raw sequencing reads.
- Generate genome references and map reads using STARsolo.
- Create and analyze Seurat objects for clustering and visualization.

### Dependencies:
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
   cd nextflow_501Project/
   ```
   
2. Make sure Nextflow is installed.\
   [Installation guide](https://www.nextflow.io/docs/latest/install.html)

4. Make sure Singularity is installed and configured. If you are using a BCGSC computer, there is no need for installation or configuration.\
   [Installation guide](https://docs.sylabs.io/guides/3.0/user-guide/installation.html)
   
6. Download the reference genome files:\
   There is a `download` process in my Nextflow pipeline that will handle the downloads automatically.

---

### Running the Pipeline

   ```bash
   cd nextflow_501Project/
   nextflow run nextflow_version.nf 
   ```

### Parameters

| Parameter            | Description                                           | Example                                    |
|----------------------|-------------------------------------------------------|--------------------------------------------|
| `fastq_dir`          | Directory containing FASTQ files                     | `data/`                                    |
| `output_dir`         | Directory for storing output files                   | `results/`                                 |
| `genome_dir`         | Directory for storing genome files                   | `ref/STAR_genomeDir/`                      |
| `threads`            | Number of threads for parallel execution             | `8`                                        |

---

## Input

- **FASTQ files**: Paired-end sequencing reads named as `*_1.fastq.gz` and `*_2.fastq.gz`.
- **Genome files**:
  - Fasta files: `chr1.fa`, `chr2.fa`, `chr3.fa`.
  - Annotation file: `gencode.v38.annotation.gtf`.

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

--- 

This README follows a structured approach to guide users in understanding and running your pipeline effectively. If you have further customizations, feel free to integrate them.
