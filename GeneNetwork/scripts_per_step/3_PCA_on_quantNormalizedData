#!/bin/bash
#
#SBATCH --job-name=run_step_3
#SBATCH --output=run_step_3.log
#SBATCH --error=run_step_3.log
#
#SBATCH --ntasks=1
#SBATCH --time=6-23:59:00
#SBATCH --mem-per-cpu=100000

module load Java/1.8.0_144-unlimited_JCE

java -jar -Xms99G -Xmx99G jars/RunV12.jar ~/brain_eQTL/GeneNetwork/steps_config_files/3_PCA_on_quantNormalizedData.sh
