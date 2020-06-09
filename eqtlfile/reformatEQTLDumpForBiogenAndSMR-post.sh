#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=6:00:00
#SBATCH --mem=15g
#SBATCH --cpus-per-task=1

ml Java
ml Python


if [ "$#" -ne 4 ]; then
    echo "Usage: dumpdir outprefix fdrreference-FDR0.05.txt.gz chr"
    exit
fi

dumpdir=$1
outprefix=$2
fdrref=$3
chr=$4

	# reformat to biogen preferred format
#	python /groups/umcg-biogen/tmp03/tools/biogen/reformatQTLfile.py \
#		$dumpdir/splitByChr/$outprefix\-$chr.txt.gz \
#		$dumpdir/biogenformat/$outprefix\-$chr\-biogenformat.txt.gz

	# grep significant hits
	# python /groups/umcg-biogen/tmp03/tools/brain_eQTL/eqtlfile/select_fdr_significant.py \
	#	$dumpdir/biogenformat/ \
	#	$dumpdir/biogenformat/fdr005/
	
	# convert to SMR file
	java -Xmx10g -jar /groups/umcg-biogen/tmp03/tools/eqtl-mapping-pipeline-1.4.8-SNAPSHOT/eqtl-mapping-pipeline.jar --mode util --converteqtlfiletosmr \
		--in $dumpdir/splitByChr/$outprefix\-$chr.txt.gz \
		--out $dumpdir/smr/$outprefix\-$chr\-SMR.txt.gz \
		--in2 $dumpdir/SNPQCLog-MAF.txt.gz
	
	# remove liftover issues
	python /groups/umcg-biogen/tmp03/tools/brain_eQTL/eqtlfile/smrFilterLiftOverIssues.py \
		$dumpdir/smr/$outprefix\-$chr\-SMR.txt.gz \
		$dumpdir/smr/$outprefix\-$chr\-SMR-fix.txt.gz \
		$chr
	rm $dumpdir/smr/$outprefix\-$chr\-SMR.txt.gz

	# remove possible duplicate eQTLs
	python /groups/umcg-biogen/tmp03/tools/brain_eQTL/eqtlfile/smrRemoveDuplicates.py \
		$dumpdir/smr/$outprefix\-$chr\-SMR-fix.txt.gz \
		$dumpdir/smr/$outprefix\-$chr\-SMR-fix-dedup.txt.gz
	rm $dumpdir/smr/$outprefix\-$chr\-SMR-fix.txt.gz
	gunzip $dumpdir/smr/$outprefix\-$chr\-SMR-fix-dedup.txt.gz

	# convert to SMR format
	/groups/umcg-wijmenga/tmp03/users/umcg-aclaringbould/Brisbane/SMR/tools/smr_Linux \
		--qfile $dumpdir/smr/$outprefix\-$chr\-SMR-fix-dedup.txt \
		--make-besd-dense \
		--out $dumpdir/smr/$outprefix\-$chr\-SMR-besd
	rm $dumpdir/smr/$outprefix\-$chr\-SMR-fix-dedup.txt


