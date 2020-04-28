#!/bin/bash
#SBATCH --job-name=REPLACENAME
#SBATCH --output=REPLACENAME.out
#SBATCH --error=REPLACENAME.err
#SBATCH --time=23:59:59
#SBATCH --mem REPLACEMEMb
#SBATCH --nodes 1
#SBATCH --qos=REPLACEQOS

ml Java

set -e
set -u

mkdir -p $(dirname REPLACEOUT)

remove=false
if [ ! -f REPLACEEIGENVECTORS ]  && [ -f REPLACEEIGENVECTORS.gz ];
then
    echo "Eigenvectors are gzipped but this step needs unzipped eigenvectors. Unzipping now"
    zcat REPLACEEIGENVECTORS.gz > REPLACEEIGENVECTORS
    remove=true
fi

# have to make sure that the input background genes are only genes that are also in eigenvectors
comm -1 -2 <(awk '{print $1}' REPLACEEIGENVECTORS | sort) <(sort REPLACEBACKGROUND) > $(basename REPLACEBACKGROUND.onlyInEigenvector.txt)

java -jar -XmsREPLACEMEM -XmxREPLACEMEM REPLACEGENENETWORKDIR/GeneNetworkBackend-1.0.7-SNAPSHOT-jar-with-dependencies.jar \
  -e REPLACEEIGENVECTORS \
  -p REPLACEIDENTITYMATRIX \
  -b $(basename REPLACEBACKGROUND.onlyInEigenvector.txt) \
  -o REPLACEOUT

touch REPLACENAME.finished



grep AUC REPLACENAME.out | cut -f2,7,11,13 >> $(dirname REPLACEOUT)/REPLACENAME.AUC.txt;
