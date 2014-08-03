default: bedtools-2.17.0/bin/bedtools boost_1_45_0_subset/boost samtools-0.1.19/samtools

bedtools-2.17.0/bin/bedtools:
	echo "preparing bedtools.v2.17.0 ..."
	tar zxvf BEDTools.v2.17.0.tar.gz
	cd bedtools-2.17.0; make
	cd ..

boost_1_45_0_subset/boost:
	echo "preparing required boost headers ..."
	tar zxvf boost_1_45_0_subset.tar.gz

samtools-0.1.19/samtools:
	echo "preparing samtools v0.1.19 ..."
	tar zxvf samtools-0.1.19.tar.gz
	cd samtools-0.1.19; make
	cd ..
