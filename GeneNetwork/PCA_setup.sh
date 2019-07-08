#!/bin/bash

source /groups/umcg-wijmenga/tmp04/projects/2019-comethylationnetwork/tools/pcapp/mkl/compilers_and_libraries_2019.3.199/linux/mkl/bin/mklvars.sh intel64
export DYLD_LIBRARY_PATH=/groups/umcg-wijmenga/tmp04/projects/2019-comethylationnetwork/tools/pcapp/mkl/compilers_and_libraries_2019.3.199/linux/mkl/lib/intel64/:/groups/umcg-wijmenga/tmp04/projects/2019-comethylationnetwork/tools/pcapp/clang/llvm-project/build/lib/:$DYLD_LIBRARY_PATH
export LD_LIBRARY_PATH=/groups/umcg-wijmenga/tmp04/projects/2019-comethylationnetwork/tools/pcapp/clang/llvm-project/build/lib/:$LD_LIBRARY_PATH

# After, can run PCA with
# /groups/umcg-wijmenga/tmp04/projects/2019-comethylationnetwork/tools/pcapp/pcaplusplus/pca.sh
