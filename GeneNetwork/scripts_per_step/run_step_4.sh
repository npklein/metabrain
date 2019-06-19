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
module load RPlus

set -e
set -u

echo "calculate geoMean"
Rscript ~/brain_eQTL/GeneNetwork/misc/calculate_geoMean.R \
    -e /groups/umcg-biogen/tmp04/umcg-ndeklein/GeneNetwork/output/samples500000reads/2019-06-12.all-datasets.FIXEDSAMPLEHEADER.kallistoAbunds_extractedColumns.txt.gz \
    -o /groups/umcg-biogen/tmp04/umcg-ndeklein/GeneNetwork/output/step3_step4/geoMean.txt

echo "Done, run step 4"
java -jar -Xms99G -Xmx99G jars/RunV12.jar ~/brain_eQTL/GeneNetwork/steps_config_files/4_DeseqNormalizedData.sh
