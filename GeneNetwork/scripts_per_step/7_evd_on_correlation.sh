#!/bin/bash
#SBATCH --job-name=evd_on_correlation
#SBATCH --output=7_evd_on_correlation.out
#SBATCH --error=7_evd_on_correlation.err
#SBATCH --time=23:59:59
#SBATCH --cpus-per-task REPLACETHREADS
#SBATCH --mem REPLACEMEMb
#SBATCH --nodes 1
#SBATCH --qos=REPLACEQOS

set -e
set -u
ml Python/3.6.3-foss-2015b
ml R
n_PCs=100

output_dir=REPLACEOUTDIR
mkdir -p $output_dir/7_evd_on_correlation_matrix/

# step 7. Run evd on correlation matrix
cd $output_dir/7_evd_on_correlation_matrix

python REPLACESCRIPTDIR/evd.py REPLACECORMATRIX REPLACECOVCORRECTEDEXPRESSION ./ --svd_solver REPLACESVDSOLVER

R --vanilla << "EOF"
eigen <- read.table('eigenvalues.txt', header=T)
png("eigenvalues.png")
plot(eigen$eigenvalues, type = "b", pch = 19, ylab = "eigenvalues")
dev.off()
EOF

mv eigenvectors.txt $output_dir/7_evd_on_correlation_matrix/REPLACENAME.eigenvectors.txt
mv pc-scores.txt $output_dir/7_evd_on_correlation_matrix/REPLACENAME.pc-scores.txt
gzip $output_dir/7_evd_on_correlation_matrix/REPLACENAME.eigenvectors.txt
gzip $output_dir/7_evd_on_correlation_matrix/REPLACENAME.pc-scores.txt

