#!/bin/bash

ml Java
ml R


coloctools=/groups/umcg-biogen/tmp03/tools/COLOCTools.jar # java coloc tools for merging eqtl files, calculating beta's etc
colocscript=/groups/umcg-biogen/tmp03/tools/brain_eQTL/colocalization/COLOC/runcoloc.r			# R script for performing coloc
eqtlfile=/groups/umcg-biogen/tmp03/output/2019-02-25-FreezeTwo/output/2019-07-08-Cortex-Cis/EUR-run2-dump/eQTLDump-sort-FDR-stripped-significantgenes.txt.gz				# eQTLDump-sort-FDR-stripped.txt.gz or eQTLDump-sort-FDR-stripped-significantgenes.txt.gz
eqtlsummarystats=/groups/umcg-biogen/tmp03/output/2019-02-25-FreezeTwo/output/2019-07-08-Cortex-Cis/EUR-run2-dump/2019-07-08-Cortex-Cis-SNPQCLog-MAF-sort.txt.gz		# 2019-07-08-Cortex-Cis-SNPQCLog-MAF-sort.txt.gz
eqtlsamples=3781		# nr of samples in the eQTL study (considered a constant, I know...)
gwasinput=$1			# PLINK formatted GWAS 
cases=$2				# nr of cases in the GWAS
controls=$3				# nr of controls in the GWAS
colocdir=$4				# outputdir

mkdir -p $colocdir/tmp/

# merge GWAS with eQTL file, output a file per gene
java -Xmx10g -jar $coloctools \
	--mergegwaswitheqtls \
	--eqtl $eqtlfile \
	--eqtlsummarystats $eqtlsummarystats \
	--gwas $gwasinput
	--output $colocdir/tmp/

# run coloc
Rscript --vanilla $colocscript $colocdir/tmp/ $eqtlsamples $cases $controls

# merge COLOC output
java -Xmx10g -jar $coloctools \
	--colocparser \
	--input $colocdir/tmp/ \
	--colocsummaryout $colocdir/$gwasinput\-coloc-summary.txt.gz \
	--colocgeneoutput $colocdir/$gwasinput\-coloc-genes.txt.gz 
	
# cleanup
rm -rv $colocdir/tmp/