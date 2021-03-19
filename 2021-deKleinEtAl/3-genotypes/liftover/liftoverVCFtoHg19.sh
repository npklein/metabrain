#ml VCFtools
#ml HTSlib


vcfin=$1
bedloc=$2
chr=$3
vcfout=$4

# convert VCF to bed
java -Xmx9g -jar /groups/umcg-biogen/tmp03/tools/VCFUtils.jar --vcftobed -i $vcfin -o $bedloc

# lift over
/groups/umcg-biogen/tmp03/tools/liftover/liftOver \
        $bedloc \
        /groups/umcg-biogen/tmp03/tools/liftover/hg38ToHg19.over.chain.gz \
        $bedloc\_hg19_lifted.bed \
        $bedloc\_hg19_unlifted.bed

		
# replace coordinates in VCF
java -Xmx9g -jar /groups/umcg-biogen/tmp03/tools/VCFUtils.jar --liftover \
		-i $vcfin \
		-b $bedloc\_hg19_lifted.bed \
		--chr $chr \
		-o $vcfout\_tmp.vcf.gz

rm $bedloc
rm $bedloc\_hg19_lifted.bed
rm $bedloc\_hg19_unlifted.bed

zcat $vcfout\_tmp.vcf.gz | vcf-sort -c | bgzip > $vcfout
tabix -p vcf $vcfout

rm $vcfout\_tmp.vcf.gz
