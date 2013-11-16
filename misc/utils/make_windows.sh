#!/bin/bash
#$ -S /bin/bash
#$ -cwd

tumor_pileup=$1
normal_pileup=$2
outf=$3
window_size=100

path=`dirname $0`

ruby $path/search_variants.rb $tumor_pileup $normal_pileup $outf/
bedtools window -v -a $outf/cand_somatic -b $outf/cand_hetero_germline -w $window_size > $outf/tmp_v_window
bedtools window -a $outf/cand_somatic -b $outf/cand_hetero_germline -w $window_size > $outf/tmp_window
cat $outf/tmp_v_window $outf/tmp_window | ruby $path/to_windows.rb | sort -k 1,1 -k 2,2n -k 3,3n > $outf/windows
rm $outf/tmp_v_window $outf/tmp_window
