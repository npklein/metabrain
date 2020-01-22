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
chainfile=/groups/umcg-biogen/tmp03/tools/liftover/hg19ToHg38.over.chain.gz

removechr=/groups/umcg-biogen/tmp03/tools/brain_eQTL/liftover/removechr.py
liftoverutils=/groups/umcg-biogen/tmp03/tools/brain_eQTL/liftover/liftOverUtils.py
ucscliftover=/groups/umcg-biogen/tmp03/tools/liftover/liftOver
allelewriter=/groups/umcg-biogen/tmp03/tools/TriTyperSNPAlleleWriter.jar

#if [ 1 -eq 2 ]; then
mkdir -p $liftoverdir
rm -v $liftoverdir/*

# strip chr label
python $removechr $indir/SNPs.txt.gz $indir/SNPs.tmp.txt.gz
python $removechr $indir/SNPMappings.txt.gz $indir/SNPMappings.tmp.txt.gz

# convert to bed, write alleles
java -Xmx8g -jar $allelewriter \
	$indir \
	$outdir/liftover/snps.bed

# lift over
$ucscliftover \
        $outdir/liftover/snps.bed \
        $chainfile \
        $outdir/liftover/snps_hg38_lifted.bed \
        $outdir/liftover/snps_hg38_unlifted.bed

# gzip
gzip -v $outdir/liftover/*.bed

# rewrite SNPMappings.txt.gz
python $liftoverutils \
	--mode updatetrityperpos \
	--bed $outdir/liftover/snps.bed.gz \
	--bedlifted $outdir/liftover/snps_hg38_lifted.bed.gz \
	--output $outdir/liftover/
#fi
# update RS ids
python $liftoverutils \
	--mode updatetrityperrsid \
	--snpmap $outdir/liftover/SNPMappings_lifted.txt.gz \
	--vcf $dbsnp \
	--output $outdir/liftover/


# liftover done
