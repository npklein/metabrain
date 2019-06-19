echo "source /home/umcg-ndeklein/brain_eQTL/GeneNetwork/PCA_setup_calculon.sh"
. /home/umcg-ndeklein/brain_eQTL/GeneNetwork/PCA_setup_calculon.sh
echo "done"
cd /groups/umcg-biogen/tmp04/umcg-ndeklein/GeneNetwork/output/step3/
echo "run PCA"
/groups/umcg-wijmenga/tmp04/projects/2019-comethylationnetwork/tools/pcapp/pcaplusplus/pca pca correlation $1
