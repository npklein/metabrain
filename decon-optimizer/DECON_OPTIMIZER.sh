#!/bin/bash
#SBATCH --job-name=DECON_OPTIMIZER
#SBATCH --output=/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-11-10-decon-optimizer/output/DECON_OPTIMIZER.out
#SBATCH --error=/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-11-10-decon-optimizer/output/DECON_OPTIMIZER.out
#SBATCH --time=23:59:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

module load Python/3.6.3-foss-2015b
source $HOME/venv/bin/activate

python3 /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-11-10-decon-optimizer/decon_optimizer.py -eq ../2020-10-12-deconvolution_gav/matrix_preparation/cortex_eur_cis/combine_eqtlprobes/eQTLprobes_combined.txt.gz -ge ../2020-10-12-deconvolution_gav/matrix_preparation/cortex_eur_cis/create_matrices/genotype_table.txt.gz -al ../2020-10-12-deconvolution_gav/matrix_preparation/cortex_eur_cis/create_matrices/genotype_alleles.txt.gz -ex ../2020-10-12-deconvolution_gav/matrix_preparation/cortex_eur_cis/create_matrices/expression_table.txt.gz -cf ../2020-10-12-deconvolution_gav/matrix_preparation/cortex_eur_cis/perform_deconvolution/deconvolution_table.txt -de ../2020-10-12-deconvolution_gav/2020-11-20-decon-QTL/cis/cortex/decon_out/deconvolutionResults.csv -sa /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-02-03-phenotype-table/2020-09-04.brain.phenotypes.withReannotatedDiagnosis.txt -sid rnaseq_id -cid cohort -a 0.01 -t MEDIUM

deactivate