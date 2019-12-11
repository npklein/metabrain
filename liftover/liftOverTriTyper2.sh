#!/bin/bash
#SBATCH --job-name=liftover
#SBATCH --output=liftover.out
#SBATCH --error=liftover.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=50gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

ml Python/2.7.11-foss-2015b
ml Java

indir=$1

outdir=$indir
liftoverdir=$outdir/liftover/
dbsnp=/groups/umcg-biogen/tmp03/annotation/All_20180418.vcf.gz
dbsnp=/groups/umcg-biogen/tmp03/annotation/SNPMappings-dbsnp151-b38-20180418.vcf.gz

#if [ 1 -eq 2 ]; then
mkdir -p $liftoverdir
rm -v $liftoverdir/*

# strip chr label
python /groups/umcg-biogen/tmp03/tools/liftover/removechr.py $indir/SNPs.txt.gz $indir/SNPs.tmp.txt.gz
python /groups/umcg-biogen/tmp03/tools/liftover/removechr.py $indir/SNPMappings.txt.gz $indir/SNPMappings.tmp.txt.gz

# convert to bed, write alleles
java -Xmx8g -jar /groups/umcg-biogen/tmp03/tools/TriTyperSNPAlleleWriter.jar \
	$indir \
	$outdir/liftover/snps.bed

# lift over
/groups/umcg-biogen/tmp03/tools/liftover/liftOver \
        $outdir/liftover/snps.bed \
        /groups/umcg-biogen/tmp03/tools/liftover/hg19ToHg38.over.chain.gz \
        $outdir/liftover/snps_hg38_lifted.bed \
        $outdir/liftover/snps_hg38_unlifted.bed

# gzip
gzip -v $outdir/liftover/*.bed

# rewrite SNPMappings.txt.gz
python /groups/umcg-biogen/tmp03/tools/liftover/liftOverUtils.py \
	--mode updatetrityperpos \
	--input $outdir/liftover/snps.bed.gz \
	--input2 $outdir/liftover/snps_hg38_lifted.bed.gz \
	--output $outdir/liftover/
#fi
# update RS ids
python /groups/umcg-biogen/tmp03/tools/liftover/liftOverUtils.py \
	--mode updatetrityperrsid \
	--input $outdir/liftover/SNPMappings_lifted.txt.gz \
	--input2 $dbsnp \
	--output $outdir/liftover/


# liftover done
