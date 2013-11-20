#!/bin/sh

tumor_pileup=$1
normal_pileup=$2
outf=$3
window_size=100

path=$(cd $(dirname $0); pwd)
bedtools_=$path/../libs/bedtools-2.17.0/bin/bedtools
mkdir -p $outf

echo "tumor pileup file: $1"
echo "normal pileup file: $2"
echo "output folder: $3"
echo "-------------------------"
echo "searching candidate mutations and hetero germline variants..."
ruby $path/search_variants.rb $tumor_pileup $normal_pileup $outf/
echo "done."
echo "search hetero germline varaints nearby for each candidate mutation..."
if [ -s $outf/cand_somatic ]; then
  if [ -s $outf/cand_hetero_germline ]; then
    $bedtools_ window -v -a $outf/cand_somatic -b $outf/cand_hetero_germline -w $window_size > $outf/tmp_v_window
    $bedtools_ window -a $outf/cand_somatic -b $outf/cand_hetero_germline -w $window_size > $outf/tmp_window
    cat $outf/tmp_v_window $outf/tmp_window | ruby $path/to_windows.rb | sort -k 1,1 -k 2,2n -k 3,3n > $outf/windows
    rm $outf/tmp_v_window $outf/tmp_window
  else
    echo "found no hetero germline variants"
    cat $outf/cand_somatic | ruby $path/to_windows.rb | sort -k 1,1 -k 2,2n -k 3,3n > $outf/windows
  fi
  echo "done.\ngenerated a list of windows: ${outf}/windows"
else
  echo "found no candidate mutations. could not make any windows."
fi
