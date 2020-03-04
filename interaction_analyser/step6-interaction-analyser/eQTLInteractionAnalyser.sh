#!/bin/bash
#SBATCH --job-name=eQTLInteractionAnalyser
#SBATCH --output=/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-03-interaction-analyser/step6-interaction-analyser/eQTLInteractionAnalyser.out
#SBATCH --time=05:59:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

module load Java

java -jar /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-03-interaction-analyser/step6-interaction-analyser/eQTLInteractionAnalyser-1.2-SNAPSHOT-jar-with-dependencies.jar \
    -input /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-03-interaction-analyser/step5-prepare-ia-inputs/output/binary/ \
    -output /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-03-interaction-analyser/step6-interaction-analyser/output/ \
    -eqtls /groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-02-18-eqtls/cortex-cis-EURandAFR-iterative/Iteration1/eQTLProbesFDR0.05-ProbeLevel.txt.gz \
    -maxcov 1 \
    -noNormalization \
    -nnoCovNormalization