TEMP_DIR=./data/temp

mkdir $TEMP_DIR
minimap2 -a ./data/fasta/dd2.fasta ./data/fastq/G22_Control_C_DD2_DRAG1_PfMAP.fastq > $TEMP_DIR/alignment.sam
samtools view -@ 4 -Sb -o $TEMP_DIR/example_alignment.bam $TEMP_DIR/alignment.sam
samtools sort -O bam -o $TEMP_DIR/sorted_example_alignment.bam $TEMP_DIR/example_alignment.bam
samtools index $TEMP_DIR/sorted_example_alignment.bam
samtools view $TEMP_DIR/sorted_example_alignment.bam > ./data/output/sorted_example_alignment.txt

rm -rf $TEMP_DIR