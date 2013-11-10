#!/bin/sh
i=$1
echo "--- test #$i ---"
cd ~/thesis/hapmuc3
rm -r hapmuc_testcases/$i/result
mkdir hapmuc_testcases/$i/result
#time ./hapmuc --tumorBamFile hapmuc_testcases/$i/tumor.sort.bam --normalBamFile hapmuc_testcases/$i/normal.sort.bam --ref /Users/usuyama/data_store/hg19_bwa-0.5.10/hg19.fasta --varFile hapmuc_testcases/$i/windows --libFile hapmuc_testcases/$i/cigar.libraries.txt --outputFile hapmuc_testcases/$i/result/mc --mapQualThreshold 0.01 > hapmuc_testcases/$i/result/log
time ./hapmuc --tumorBamFile hapmuc_testcases/$i/tumor.sort.bam --normalBamFile hapmuc_testcases/$i/normal.sort.bam --ref /Users/usuyama/data_store/hg19_bwa-0.5.10/hg19.fasta --varFile hapmuc_testcases/$i/windows --outputFile hapmuc_testcases/$i/result/mc --mapQualThreshold 0.9683772233983162 > hapmuc_testcases/$i/result/log
cat hapmuc_testcases/$i/README
tail hapmuc_testcases/$i/windows
tail hapmuc_testcases/$i/windows | ruby -ane 'puts $F[-2..-1].join("\t")'
tail -n1 hapmuc_testcases/$i/result/mc.calls.txt | ruby -ane 'puts $F[-3..-1].join("\t")'
