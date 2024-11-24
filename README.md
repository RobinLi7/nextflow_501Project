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

---

## Usage

### Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/RobinLi7/nextflow_501Project.git
   cd nextflow_501Project/
   ```

2. Install **Nextflow**:
   ```bash
   curl -s https://get.nextflow.io | bash
   ```

3. Make sure Singularity is installed and configured:
   ```bash
   sudo apt-get install singularity
   ```

4. Download the reference genome files:
   Use the `download` process or manually run:
   ```bash
   curl -L -o ref/CellRanger/3M-february-2018.txt.gz https://github.com/noamteyssier/10x_whitelist_mirror/raw/main/3M-february-2018.txt.gz
   curl -L -o ref/chr1.fa.gz http://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr1.fa.gz
   ```
   (Additional download commands included in the pipeline.)

---

### Running the Pipeline

1. Prepare input FASTQ files:
   Organize paired-end FASTQ files in `data/` directory. Each pair should follow the naming convention `*_1.fastq.gz` and `*_2.fastq.gz`.

2. Configure parameters:
   Edit `nextflow.config` to specify the directories for input and output.

3. Execute the pipeline:
   ```bash
   nextflow run main.nf -params-file params.yaml
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
