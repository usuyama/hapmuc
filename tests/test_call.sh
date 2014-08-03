#!/bin/bash
#$ -S /bin/bash
#$ -cwd

i=$1
echo "index: $1"
mkdir -p $i/result
../bin/hapmuc -a data/$i/tumor.bam -b data/$i/normal.bam -f data/random_ref.fasta -w data/$i/windows -o data/$i/result/mc
echo "please check data/$i/result/mc.calls.txt for result"
cat data/$i/result/mc.calls.txt
