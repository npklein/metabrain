#!/bin/bash
#SBATCH --job-name=PCA
#SBATCH --output=PCA.out
#SBATCH --error=PCA.err
#SBATCH --time=23:59:59
#SBATCH --cpus-per-task 22
#SBATCH --mem 100gb
#SBATCH --nodes 1
#SBATCH --qos=dev

source REPLACEPCASETUP
output_dir=REPLACEOUTDIR
outfile_step6=REPLACEOUTFILE
mkdir -p $output_dir/7_PCA_on_correlation_matrix/
if [ ! -f $output_dir/7_PCA_on_correlation_matrix/pc-scores.txt ];
then
    # step 7. Run PCA on correlation matrix
#    cd $TMPDIR
    cd $output_dir/7_PCA_on_correlation_matrix
    /groups/umcg-wijmenga/tmp04/projects/2019-comethylationnetwork/tools/pcapp/pcaplusplus/pca evd $outfile_step6
    /groups/umcg-wijmenga/tmp04/projects/2019-comethylationnetwork/tools/pcapp/pcaplusplus/pca pc-scores correlation $outfile_step6 eigenvectors.txt

#    mv summary.txt $output_dir/7_PCA_on_correlation_matrix/MetaBrain.summary.txt
#    mv correlation.txt $output_dir/7_PCA_on_correlation_matrix/MetaBrain.correlation.txt
    mv eigenvalues.txt $output_dir/7_PCA_on_correlation_matrix/MetaBrain.eigenvalues.txt
    mv eigenvectors.txt $output_dir/7_PCA_on_correlation_matrix/MetaBrain.eigenvectors.txt
    mv pc-scores.txt $output_dir/7_PCA_on_correlation_matrix/MetaBrain.pc-scores.txt
    mv cronbach.txt $output_dir/7_PCA_on_correlation_matrix/MetaBrain.cronbach.txt.txt
    cd -
fi

