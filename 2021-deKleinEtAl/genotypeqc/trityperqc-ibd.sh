ml Java
ml plink
ml Python/2.7.12-foss-2015b

dsname=$1
dsdir=$2

snpmapref=/groups/umcg-biogen/tmp03/annotation/GiantSNPMappings_Filtered_on_Maf0.001_Indels.txt.gz
kgbim=/groups/umcg-biogen/tmp03/annotation/1kgp3v5a/ALL-merged.bim
kgpedmap=/groups/umcg-biogen/tmp03/annotation/1kgp3v5a/ALL-merged


vcfutils=/groups/umcg-biogen/tmp03/tools/VCFUtils.jar
ttvcfconverter=/groups/umcg-biogen/tmp03/tools/ConvertTTToVCF.jar
updatersscript=/groups/umcg-biogen/tmp03/tools/genotypeqc/updaters.py
filteranddedupscript=/groups/umcg-biogen/tmp03/tools/liftover/filteranddedup.py
missingnessplotscript=/groups/umcg-biogen/tmp03/tools/lmiss-hist.Rscript

workdir=$dsdir/qcout/
mkdir -p $workdir

# convert to VCF
java -Xmx4g -jar $ttvcfconverter \
	$dsdir \
	$workdir \
	0.05 0.9 0.0001

#read -p "Press enter to continue1"
# convert to plink
plink --vcf $workdir/ConvertedFromTriTyper.vcf.gz \
	--recode \
	--const-fid 0 \
	--make-bed --out $workdir/$dsname

plink --bfile $workdir/$dsname \
	--const-fid \
	--indep-pairwise 25000 2000 0.2 --out $workdir/$dsname\-fulldata

# MDS
plink --bfile $workdir/$dsname --genome \
	--const-fid \
	--extract $workdir/$dsname\-fulldata.prune.in \
	--out $workdir/$dsname\-genome-ibd

plink --bfile $workdir/$dsname --read-genome $workdir/$dsname\-genome-ibd.genome \
	--const-fid \
	--extract $workdir/$dsname\-fulldata.prune.in \
	--cluster --mds-plot 4 --out $workdir/$dsname\-mds-ibd

rm $workdir/*.bed
rm $workdir/*.bim
rm $workdir/*.fam
rm $workdir/*.vcf.gz
rm $workdir/*.ped
rm $workdir/*.map
rm $workdir/ConversionLog.txt.gz
rm $workdir/*.prune.in
rm $workdir/*.prune.out
rm $workdir/*.nosex
rm $workdir/*.log

gzip -v $workdir/*
