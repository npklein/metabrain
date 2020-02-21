#!/bin/bash
module load Python/2.7.12-foss-2015b
module load Java/11.0.2

indir=$1
output=$2

summarystats=/groups/umcg-wijmenga/tmp03/projects/eQTLGen/helpfiles/ePRSHelpFiles/SummaryStatisticsDatabase_freeze1b/

mkdir -p $output

# fix snp ids
mkdir -p $output/tmpgt/

cp -v $indir/* $output/tmpgt/

# kill b38 data
python /groups/umcg-biogen/tmp03/tools/prstools/stripAndDedup.py $output/tmpgt/SNPs.txt.gz

# move in hg19 snp mappings
cp -v /groups/umcg-wijmenga/tmp03/projects/eQTLGen/helpfiles/GiantSNPMappings_Filtered_on_Maf0.001_Indels.txt.gz $output/tmpgt/SNPMappings.txt.gz

java -Xmx80g -Xms80g -jar /groups/umcg-biogen/tmp03/tools/GeneticRiskScoreCalculator-0.1.4.1-SNAPSHOT/GeneticRiskScoreCalculator.jar \
	-gi $output/tmpgt/ \
	-gt TRITYPER \
	-i $summarystats \
	-o $output \
	-p 5e-8:1e-5:1e-4:1e-3:1e-2 -r 0.1 \
	-w 250000:5000000 \
	-er 6:25000000-35000000 2>&1 | tee $output/PRSCalculation.log

#rm -rv $output/tmpgt/
#gzip -v $output/*.txt

