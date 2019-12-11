#!/bin/bash
module load Python/2.7.12-foss-2015b
module load Java/11.0.2
module load plink/1.9-foss-2015b


indir=$1
output=$2


mkdir -p $output

# fix snp ids
mkdir -p $output/tmpgt/

cp -v $indir/* $output/tmpgt/

# kill b38 data
python /groups/umcg-biogen/tmp03/tools/prstools/stripAndDedup.py $output/tmpgt/SNPs.txt.gz

# move in hg19 snp mappings
cp -v /groups/umcg-wijmenga/tmp03/projects/eQTLGen/helpfiles/GiantSNPMappings_Filtered_on_Maf0.001_Indels.txt.gz $output/tmpgt/SNPMappings.txt.gz

# convert to vcf
mkdir $output/vcf/
java -Xmx20g -jar /groups/umcg-biogen/tmp03/tools/ConvertTTToVCF.jar \
	$output/tmpgt/ \
	$output/vcf/ \
	0 0.5 0.0001

#rm -rv $output/tmpgt/
mkdir $output/plink/

# convert to plink
plink --vcf $output/vcf/ConvertedFromTriTyper.vcf.gz \
	--const-fid 0 \
	--make-bed --out $output/plink/dataset

#rm -rv $output/vcf/
mkdir $output/plinkfix/

# make sure chromosome problems are solved
plink --bfile $output/plink/dataset \
	--const-fid 0 \
	--make-bed --out $output/plinkfix/dataset

#rm -rv $output/plink/

# calculate PRSs
summarystats=/groups/umcg-biogen/tmp03/output/2019-02-25-FreezeTwo/input/2019-06-18-PRS/ukbb_prs/beta_prs/*.prsgwassumstat
for stat in $summarystats;
do

echo $stat
echo $(basename "$stat".txt)
bs=$(basename "$stat".txt)

/groups/umcg-biogen/tmp03/tools/plink2 --allow-no-sex \
        --bfile $output/plinkfix/dataset \
        --score $stat ignore-dup-ids no-mean-imputation cols=nmissallele,dosagesum,scoreavgs,scoresums \
        --out $output/$bs --threads 16
done
gzip $output/*

