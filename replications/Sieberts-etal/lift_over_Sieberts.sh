echo "make bedfile"


echo "start liftover"
/groups/umcg-biogen/tmp03/tools/liftover/liftOver \
       data/Sieberts.uniqueSNPs.bed \
       /groups/umcg-biogen/tmp03/tools/liftover/hg19ToHg38.over.chain.gz \
       data/Sieberts_lifted.bed \
       data/Sieberts_unlifted.bed

echo "Done"
