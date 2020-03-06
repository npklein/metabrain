#!/bin/bash
#SBATCH --job-name=eQTLInteractionAnalyserPerGroup
#SBATCH --output=/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-03-interaction-analyser/step6-interaction-analyser/eQTLInteractionAnalyserPerGroup.out
#SBATCH --time=05:59:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

module load Java
module load Python/3.6.3-foss-2015b

/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-03-interaction-analyser/step6-interaction-analyser/eqtl_interaction_analyser_per_group.py
