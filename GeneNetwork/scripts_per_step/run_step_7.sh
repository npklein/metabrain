#!/bin/bash
#
#SBATCH --job-name=run_step_3
#SBATCH --output=run_step_3.log
#SBATCH --error=run_step_3.log
#
#SBATCH --ntasks=1
#SBATCH --time=6-23:59:00
#SBATCH --mem-per-cpu=100000


set -e
set -u

mkdir -p output/step7/
/home/umcg-ndeklein/brain_eQTL/Genash eNetwork/steps_config_files/7_PCA_on_correlationMatrix.sh

mv summary.txt output/step7/
mv covariance.txt output/step7/
mv eigenvalues.txt output/step7/
mv /eigenvectors.txt output/step7/
