#!/bin/bash
#SBATCH --job-name=evd_on_correlation
#SBATCH --output=7_evd_on_correlation.out
#SBATCH --error=7_evd_on_correlationPCA.err
#SBATCH --time=23:59:59
#SBATCH --cpus-per-task REPLACETHREADS
#SBATCH --mem 100gb
#SBATCH --nodes 1
#SBATCH --qos=regular

source REPLACEPCASETUP
set -e
set -u
output_dir=REPLACEOUTDIR
mkdir -p $output_dir/7_evd_on_correlation_matrix/

# step 7. Run evd
 on correlation matrix
cd $output_dir/7_evd_on_correlation_matrix
python REPLACESCRIPTDIR/evd.py REPLACECORMATRIX REPLACECOVCORRECTEDEXPRESSION ./ --svd_solver REPLACESVDSOLVER

n=0
while read line;
do
    if [[ $line == "Cronbach"* ]]; then continue; fi
    compare=$(echo $line'<'REPLACECRONBACHCUTOFF | bc)
    if [[ compare -eq 1 ]]; then break; fi
        ((n=$n+1))
done < cronbach.txt

cat < eigenvectors.txt | cut -f1-$n > MetaBrain.eigenvectors.cronbach_REPLACECRONBACHCUTOFF.txt
mv eigenvectors.txt $output_dir/7_PCA_on_correlation_matrix/MetaBrain.eigenvectors.txt
mv pc-scores.txt $output_dir/7_PCA_on_correlation_matrix/MetaBrain.pc-scores.txt
gzip $output_dir/7_PCA_on_correlation_matrix/MetaBrain.eigenvectors.txt
gzip $output_dir/7_PCA_on_correlation_matrix/MetaBrain.pc-scores.txt
# need unzippd file for next step, so only use zcat
zcat $output_dir/7_PCA_on_correlation_matrix/MetaBrain.eigenvectors.cronbach_REPLACECRONBACHCUTOFF.txt > $output_dir/7_PCA_on_correlation_matrix/MetaBrain.eigenvectors.cronbach_REPLACECRONBACHCUTOFF.txt.gz
cd -

