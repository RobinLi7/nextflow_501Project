# Define paths
GENOME_DIR=/projects/rli_prj/CPSC501/project/ref/STAR_genomeDir
FASTA_FILE="/projects/rli_prj/CPSC501/project/ref/chr1.fa \
               /projects/rli_prj/CPSC501/project/ref/chr2.fa \
               /projects/rli_prj/CPSC501/project/ref/chr3.fa"
GTF_FILE=/projects/rli_prj/CPSC501/project/ref/chr1_chr2_chr3.gtf

# Create genome index
STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir $GENOME_DIR \
     --genomeFastaFiles $FASTA_FILE \
     --sjdbGTFfile $GTF_FILE \
     --sjdbOverhang 90 \
     --genomeSAindexNbases 13