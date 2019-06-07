echo "make bedfile"

bedfile=data/Sieberts.bed
if [ ! -f $bedfile ];
then
    while read line;
    do
        echo "$line" | awk -F"," '{ print "chr" $1 "\t" $2 "\t" $2}' >> $bedfile
    done < 2019-06-07-Sieberts-eQTLs/Cortex_MetaAnalysis_ROSMAP_CMC_HBCC_Mayo_cis_eQTL_release.csv
fi


echo "start liftover"
/groups/umcg-biogen/tmp03/tools/liftover/liftOver \
       Sieberts.bed \
       /groups/umcg-biogen/tmp03/tools/liftover/hg19ToHg38.over.chain.gz \
       data/Sieberts_lifted.bed \
       data/Sieberts_unlifted.bed

echo "Done"
