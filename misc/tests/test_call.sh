#!/bin/bash
#$ -S /bin/bash
#$ -cwd

i=$1
echo "index: $1"
mkdir -p $i/result
../../hapmuc --tumorBamFile $i/tumor.bam --normalBamFile $i/normal.bam --ref random_ref.fasta --varFile $i/windows --outputFile $i/result/mc > $i/result/log
echo "please check $i/result/mc.calls.txt for result"
