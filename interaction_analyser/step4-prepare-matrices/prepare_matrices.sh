#!/bin/bash
#SBATCH --job-name=prepare_matrices
#SBATCH --output=/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-03-interaction-analyser/step4-prepare-matrices/prepare_matrices.out
#SBATCH --time=01:59:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

module load Python/3.4.1-foss-2015b
module load pandas/0.19.1-foss-2015b-Python-3.4.1

/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-03-interaction-analyser/step4-prepare-matrices/prepare_matrices.py
