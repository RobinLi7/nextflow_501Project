# Define paths
GENOME_DIR=/projects/rli_prj/CPSC501/project/ref/STAR_genomeDir
FASTA_FILE=/projects/rli_prj/CPSC501/project/ref/hg38.fa
GTF_FILE=/projects/rli_prj/CPSC501/project/ref/gencode.v38.annotation.gtf

# Create genome index
STAR --runThreadN 4 \
     --runMode genomeGenerate \
     --genomeDir $GENOME_DIR \
     --genomeFastaFiles $FASTA_FILE \
     --sjdbGTFfile $GTF_FILE \
     --sjdbOverhang 90
