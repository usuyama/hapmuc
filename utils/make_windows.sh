#!/bin/sh

tumor_pileup=$1
normal_pileup=$2
outf=$3
window_size=${4-600} #default=600

script_path=$(cd $(dirname $0); pwd)
bedtools_=$script_path/../libs/bedtools-2.17.0/bin/bedtools
mkdir -p $outf

echo "tumor pileup file: $1"
echo "normal pileup file: $2"
echo "output folder: $3"
echo "-------------------------"
echo "searching candidate mutations and hetero germline variants..."
python $script_path/search_variants/search_variants.py $tumor_pileup $normal_pileup $outf/

cand_somatic_file=$outf/somatic_candidates
cand_hetero_germline_file=$outf/hetero_germline_candidates

echo "done."
echo "search hetero germline varaints nearby for each candidate mutation..."
if [ -s $cand_somatic_file ]; then
  if [ -s $cand_hetero_germline_file ]; then
    $bedtools_ window -v -a $cand_somatic_file -b $cand_hetero_germline_file -w $window_size > $outf/tmp_v_window
    $bedtools_ window -a $cand_somatic_file -b $cand_hetero_germline_file -w $window_size > $outf/tmp_window
    cat $outf/tmp_v_window $outf/tmp_window | ruby $script_path/to_windows.rb | sort -k 1,1 -k 2,2n -k 3,3n > $outf/windows
    #rm $outf/tmp_v_window $outf/tmp_window
  else
    echo "found no hetero germline variants"
    cat $cand_somatic_file | ruby $script_path/to_windows.rb | sort -k 1,1 -k 2,2n -k 3,3n > $outf/windows
  fi
  echo "done.\ngenerated a list of windows: ${outf}/windows"
else
  echo "found no candidate mutations. could not make any windows."
fi
