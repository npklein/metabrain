#!/bin/bash
#SBATCH --job-name=REPLACENAME
#SBATCH --output=REPLACENAME.out
#SBATCH --error=REPLACENAME.err
#SBATCH --time=5:59:59
#SBATCH --mem 10gb
#SBATCH --nodes 1
#SBATCH --qos=regular

ml Java

set -e
set -u
mkdir -p $(dirname REPLACEOUT)
java -jar -Xms30g -Xmx30g REPLACEGENENETWORKDIR/GeneNetworkBackend-1.0.7-SNAPSHOT-jar-with-dependencies.jar \
  -e REPLACEEIGENVECTORS \
  -p REPLACEIDENTITYMATRIX \
  -b REPLACEGENENETWORKDIR/PathwayMatrix/REPLACETYPE_genesInPathways_filteredOnEigenVectorGenes.txt \
  -o REPLACEOUT

touch REPLACENAME.finished
