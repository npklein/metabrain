#!/bin/bash
#SBATCH --job-name=create_regression_matrix
#SBATCH --output=/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-03-zscore-comparison/step1-regression-matrix/create_regression_matrix.out
#SBATCH --time=00:59:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

module load Python/3.4.1-foss-2015b
module load pandas/0.19.1-foss-2015b-Python-3.4.1

/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-03-zscore-comparison/step1-regression-matrix/create_regression_matrix.py
