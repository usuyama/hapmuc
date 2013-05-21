#!/bin/sh

range=`samtools view $1.sort.bam | head -n1 | ruby -ane 'puts "#{$F[2]}:#{$F[3].to_i}-#{$F[3].to_i + 800}"'`
wd=`pwd`
base=`pwd | ruby -ane 'puts $F[0].split("/")[-3..-1].join("_")'`
cat << EOT > igv.batch
snapshotDirectory ~/tmp/igvss
new
load $wd/tumor.sort.bam
load $wd/normal.sort.bam
genome hg19
goto $range
collapse
snapshot ${base}_${range}.png
EOT
