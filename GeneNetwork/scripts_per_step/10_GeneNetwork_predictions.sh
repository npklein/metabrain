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

java -jar -XmsREPLACEMEM -XmxREPLACEMEM REPLACEGENENETWORKDIR/GeneNetworkBackend-1.0.7-SNAPSHOT-jar-with-dependencies.jar \
  -e REPLACEEIGENVECTORS \
  -p REPLACEIDENTITYMATRIX \
  -b REPLACEBACKGROUND \
  -o REPLACEOUT

touch REPLACENAME.finished


if [ "$remove" = true ];
then
    echo "Removing unzipped eigenvector file"
    rm REPLACEEIGENVECTORS
fi

