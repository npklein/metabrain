dsname=$1
mkdir -p ./qcout/

# simple maf filter
plink --bfile $dsname --mind 0.1 --maf 0.01 \
	--geno 0.05 \
	--hwe 0.00001 \
	--make-bed \
	--out ./qcout/$dsname\-filterstep1

#read -p "Press enter to continue1"

# missingness
plink --bfile $dsname \
	--missing \
	--out ./qcout/$dsname

#read -p "Press enter to continue2"

# R script for missingness
Rscript /groups/umcg-biogen/tmp03/tools/lmiss-hist.Rscript ./qcout/$dsname.lmiss ./qcout/$dsname.lmiss.pdf

#read -p "Press enter to continue3"
# heterozygosity
plink --bfile ./$dsname --het --out ./qcout/$dsname
#read -p "Press enter to continue4"
# compare to 1kg
python /groups/umcg-biogen/tmp03/tools/liftover/filteranddedup.py \
       /groups/umcg-biogen/tmp03/annotation/1kgp3v5a/ALL-merged.bim \
       ./$dsname.bim \
       ./qcout/snpintersect.txt
#read -p "Press enter to continue5"
plink --bfile $dsname --extract ./qcout/snpintersect.txt --make-bed --out ./qcout/$dsname\-1kgsnps
#read -p "Press enter to continue6"
plink --bfile /groups/umcg-biogen/tmp03/annotation/1kgp3v5a/ALL-merged \
	--extract ./qcout/snpintersect.txt --make-bed --out ./qcout/1kg-snpintersect
#read -p "Press enter to continue7"
plink --bfile ./qcout/1kg-snpintersect \
	--bmerge ./qcout/$dsname\-1kgsnps --make-bed --out ./qcout/$dsname\-1kgmerged
#read -p "Press enter to continue8"

if [ ! -f ./qcout/$dsname\-1kgmerged.bim ]; then
        plink --bfile ./qcout/$dsname\-1kgsnps --make-bed \
		--out ./qcout/$dsname\-1kgsnps-filter \
		--exclude ./qcout/$dsname\-1kgmerged-merge.missnp
#read -p "Press enter to continue9"
        plink --bfile ./qcout/1kg-snpintersect --make-bed \
		--out ./qcout/1kg-snpintersect-filter \
		--exclude ./qcout/$dsname\-1kgmerged-merge.missnp
#read -p "Press enter to continue10"
        plink --bfile ./qcout/1kg-snpintersect-filter \
		--bmerge ./qcout/$dsname\-1kgsnps-filter \
		--make-bed --out ./qcout/$dsname\-1kgmerged
#read -p "Press enter to continue11"


fi

plink --bfile ./qcout/$dsname\-1kgmerged --indep-pairwise 50 5 0.2 --out ./qcout/$dsname\-1kgmerged
#read -p "Press enter to continue12"
plink --bfile ./qcout/$dsname\-1kgmerged --extract ./qcout/$dsname\-1kgmerged.prune.in --pca 10 --out ./qcout/$dsname\-pca
#read -p "Press enter to continue13"

# IBS
plink --bfile ./qcout/$dsname\-1kgsnps --recode vcf --out ./qcout/$dsname\-1kgsnps
#read -p "Press enter to continue14"
java -Xmx20g -jar /groups/umcg-biogen/tmp03/tools/VCFUtils.jar --similarity \
        -l ./qcout/$dsname\-1kgmerged.prune.in \
        -i ./qcout/$dsname\-1kgsnps.vcf \
        -o ./qcout/$dsname\-geneticsimilarity

# MDS
plink --bfile $dsname --indep-pairwise 50 5 0.5 --out ./qcout/$dsname\-fulldata
plink --bfile $dsname --genome \
	--extract ./qcout/$dsname\-fulldata.prune.in --out ./qcout/$dsname\-genome
plink --bfile $dsname --read-genome ./qcout/$dsname\-genome.genome \
	--extract ./qcout/$dsname\-fulldata.prune.in \
	--cluster --mds-plot 4 --out ./qcout/$dsname\-mds
