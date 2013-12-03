HapMuC
======================
HapMuC is a **somatic mutation caller**, which can utilize the information of heterozygous germline variants near candidate mutations. It takes a tumor bam and a normal bam as inputs, and generates a list of candidate somatic mutations.

How to build
----------
### Prepare build dependencies ###
For convenience, the prerequisites of HapMuC are all included in this repository (You can find them at `libs/`). The prerequisites are:
* Boost 1.45.0 (modified to include only the subset of the library used by HapMuC)
* SAMtools 0.1.19
* BEDTools 2.17.0
* SeqAn (It is included as HapMuC uses its Needleman-Wunsch algorithm to align candidate haplotypes to the reference sequence.)

To prepare them, please execute the following script in your shell:
```sh
% make dependencies
```

### Build HapMuC ###

Just execute the following script in your shell:
```sh
% make
```

How to run
----------
HapMuC workflow is consisting from three steps. Let's try it on the sample bam files, which we prepared in the tests folder.
Please go to `tests/1`.
### Step1: _samtools mpileup_ ###
```sh
../../libs/samtools-0.1.19/samtools mpileup -B -f ../random_ref.fasta normal.bam > normal.pileup
../../libs/samtools-0.1.19/samtools mpileup -B -f ../random_ref.fasta tumor.bam > tumor.pileup
```
### Step2: making candidate windows ###
```sh
sh ../../utils/make_windows.sh tumor.pileup normal.pileup ./
```
Outputs:
* ./cand_somatic
    * Candidate somatic mutations, which passed the minimum criterion. (NOTE: you can set the parameters in `utils/search_variants.rb` )
* ./cand_hetero_germline
    * Heterozygous germline variants. Also, you can check the criterion in `utils/search_variants.rb`.
* ./windows
    * A list of candidate windows, which is generated based on the above two files using `bedtools window`.

### Step3: mutation calling by HapMuC algorithm ###
```sh
../../hapmuc -a tumor.bam -b normal.bam -f ../random_ref.fasta -w windows -o result/mc
```
You can check the results in `results/mc.calls.txt`

Publication
----------
Preparing for submission.

License
----------
Copyright &copy; 2013 Naoto Usuyama  
Released under the [GNU General Public License, Version 3.0][GPL].

This implementation forked from the program [Dindel][dindel], which is licensed under the GPLv3.  

[GPL]: http://www.gnu.org/licenses/gpl.html
[dindel]: http://www.sanger.ac.uk/resources/software/dindel/
