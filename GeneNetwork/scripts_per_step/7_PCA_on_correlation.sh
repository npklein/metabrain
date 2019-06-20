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
if [ ! -f $output_dir/7_PCA_on_correlation_matrix/pc-scores.txt ];
then
    # step 7. Run PCA on correlation matrix
    cd $TMPDIR
    /groups/umcg-wijmenga/tmp04/projects/2019-comethylationnetwork/tools/pcapp/pcaplusplus/pca pca correlation $outfile_step6
    /groups/umcg-wijmenga/tmp04/projects/2019-comethylationnetwork/tools/pcapp/pcaplusplus/pca pc-scores correlation $outfile_step6 eigenvectors.txt
    /groups/umcg-wijmenga/tmp04/projects/2019-comethylationnetwork/tools/pcapp/pcaplusplus/pca evd filename

    mkdir -p $output_dir/7_PCA_on_correlation_matrix/
    mv summary.txt correlation.txt eigenvalues.txt eigenvectors.txt pc-scores.txt cronbach.txt $output_dir/7_PCA_on_correlation_matrix/
    cd -
fi

