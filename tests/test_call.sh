#!/bin/bash
#$ -S /bin/bash
#$ -cwd

i=$1
echo "index: $1"
mkdir -p $i/result
../hapmuc -a $i/tumor.bam -b $i/normal.bam -f random_ref.fasta -w $i/windows -o $i/result/mc
echo "please check $i/result/mc.calls.txt for result"
cat $i/result/mc.calls.txt
