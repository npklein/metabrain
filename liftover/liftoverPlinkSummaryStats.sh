
ml Python/2.7.12-foss-2015b

# $1 : input plink file
liftoverutils=/groups/umcg-biogen/tmp03/tools/liftover/liftOverUtils.py
liftovertool=/groups/umcg-biogen/tmp03/tools/liftover/liftOver
buildout=hg38
chainfile=/groups/umcg-biogen/tmp03/tools/liftover/hg19ToHg38.over.chain.gz


python $liftoverutils --mode gwastobed \
	--in $1 --bed $1\.bed
	
$liftovertool \
	$1\.bed \
	$chainfile \
	$1\-$buildout\.bed \
	$1\-$buildout\-unlifted.bed
	
python $liftoverutils --mode liftovergwas \
	--in $1 \
	--bed $1\-$buildout.bed \
	--out $1\-$buildout.txt.gz

#rm $1\-$buildout\.bed
#rm $1\-$buildout\-unlifted.bed