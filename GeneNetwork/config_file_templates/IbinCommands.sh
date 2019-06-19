export DYLD_LIBRARY_PATH="/opt/intel/mkl/lib/":"/opt/intel/compilers_and_libraries_2016.1.111/mac/compiler/lib":$DYLD_LIBRARY_PATH
cd /groups/umcg-wijmenga/tmp03/users/umcg-svandam/skyMap2/samples500000reads/70percMapping/1_covariance/
/Volumes/Promise_RAID/juha/PCA++/pca evd /groups/umcg-wijmenga/tmp03/users/umcg-svandam/skyMap2/samples500000reads/70percMapping/1_covariance/gene_covariance.txt
mv eigenvectors.txt /groups/umcg-wijmenga/tmp03/users/umcg-svandam/skyMap2/samples500000reads/70percMapping/1_covariance/GENE.eigenvectors.txt
/Volumes/Promise_RAID/juha/PCA++/pca pc-scores covariance /groups/umcg-wijmenga/tmp03/users/umcg-svandam/skyMap2/samples500000reads/70percMapping/1_covariance/pre_Correlation_Or_Covariance.txt GENE.eigenvectors.txt
mv eigenvalues.txt /groups/umcg-wijmenga/tmp03/users/umcg-svandam/skyMap2/samples500000reads/70percMapping/1_covariance/GENE.eigenvalues.txt
n=0
while read line; do 
if [[ $line == "Cronbach"* ]]; then continue; fi
compare=$(echo $line'<'0.7 | bc)
if [[ compare -eq 1 ]]; then break; fi
((n=$n+1))
done < cronbach.txt
echo $n
cat < GENE.eigenvectors.txt | cut -f1-$n > GENE.eigenvectors0.7.txt
gzip -f GENE.eigenvectors.txt
