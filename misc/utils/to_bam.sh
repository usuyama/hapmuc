samtools view -Sb $1.sam > $1.bam
samtools sort $1.bam $1.sort
mv $1.sort.bam $1.bam
samtools index $1.bam
