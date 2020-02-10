#!/bin/bash
ml Python/2.7.12-foss-2015b
ml Java

input=$1
chainfile=/groups/umcg-biogen/tmp03/tools/liftover/hg38ToHg19.over.chain.gz
geneannotation=/groups/umcg-biogen/tmp03/annotation/Homo_sapiens.GRCh37.71-genes.txt.gz
dbsnp=/groups/umcg-biogen/tmp03/annotation/dbsnp/b37/All_20180423.vcf.gz

python /groups/umcg-biogen/tmp03/tools/liftover/liftOverUtils.py --mode eqtltobed \
	-i $input \
	-o $input\.bed


/groups/umcg-biogen/tmp03/tools/liftover/liftOver \
	$input\.bed \
	$chainfile  \
	$input\.lifted.bed \
	$input\.unlifted.bed

python /groups/umcg-biogen/tmp03/tools/liftover/liftOverUtils.py --mode liftovereqtl \
	-i $input \
	-o $input\.liftover-b37.txt.gz \
	-i2 $geneannotation \
	-i3 $input\.lifted.bed \
	-i4 $dbsnp

rm $input\.bed
rm $input\.lifted.bed
rm $input\.unlifted.bed
