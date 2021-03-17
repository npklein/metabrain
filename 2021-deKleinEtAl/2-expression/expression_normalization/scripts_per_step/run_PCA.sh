#!/bin/bash
#SBATCH --job-name=PCA
#SBATCH --output=PCA.out
#SBATCH --error=PCA.err
#SBATCH --time=23:59:59
#SBATCH --cpus-per-task 22
#SBATCH --mem 100gb
#SBATCH --nodes 1
#SBATCH --qos=dev

. ~/brain_eQTL/GeneNetwork/PCA_setup_calculon.sh
echo "done"
cd REPLACEOUT
echo "run PCA"
/groups/umcg-wijmenga/tmp04/projects/2019-comethylationnetwork/tools/pcapp/pcaplusplus/pca evd REPLACEGENECOVARIANCE
mv eigenvectors.txt REPLACEOUT/GENE.eigenvectors.txt
/groups/umcg-wijmenga/tmp04/projects/2019-comethylationnetwork/tools/pcapp/pcaplusplus/pca pc-scores covariance REPLACEGENECOVARIANCE REPLACEOUT/GENE.eigenvectors.txt
mv eigenvalues.txt REPLACEOUT/GENE.eigenvalues.txt
n=0
while read line; do
if [[ $line == "Cronbach"* ]]; then continue; fi
compare=$(echo $line'<'0.7 | bc)
if [[ compare -eq 1 ]]; then break; fi
((n=$n+1))
done < cronbach.txt
echo $n
cat < GENE.eigenvectors.txt | cut -f1-$n > GENE.eigenvectors0.7.txt
gzip -f GENE.eigenvectors.txt
