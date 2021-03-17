
echo "start liftover"
/groups/umcg-biogen/tmp03/tools/liftover/liftOver \
       $1 \
       /groups/umcg-biogen/tmp03/tools/liftover/hg19ToHg38.over.chain.gz \
       $2 \
       $3
echo "Done"
