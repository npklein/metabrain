#!/bin/bash

if [ $HOSTNAME=="boxy" ];
then
    TMP='tmp03'
elif [ $HOSTNAME=="calculon" ];
then
    TMP='tmp04'
else
    echo "HOSTNAME should be calculon or boxy, was $HOSTNAME"
    exit 1
fi
source /groups/umcg-wijmenga/$TMP/projects/2019-comethylationnetwork/tools/pcapp/mkl/compilers_and_libraries_2019.3.199/linux/mkl/bin/mklvars.sh intel64
echo "adding to DYLD_LIBRARY_PATH:"
intel64="/groups/umcg-wijmenga/$TMP/projects/2019-comethylationnetwork/tools/pcapp/mkl/compilers_and_libraries_2019.3.199/linux/mkl/lib/intel64/"
echo $intel64
llvm="/groups/umcg-wijmenga/$TMP/projects/2019-comethylationnetwork/tools/pcapp/clang/llvm-project/build/lib/"
echo $llvm

echo "adding to LD_LIBRARY_PATH:"
echo $llvm

export DYLD_LIBRARY_PATH=$intel64:$llvm:$DYLD_LIBRARY_PATH
export LD_LIBRARY_PATH=$llvm:$LD_LIBRARY_PATH

echo $DYLD_LIBRARY_PATH
echo $LD_LIBRARY_PATH
# After, can run PCA with
/groups/umcg-wijmenga/$TMP/projects/2019-comethylationnetwork/tools/pcapp/pcaplusplus/pca
