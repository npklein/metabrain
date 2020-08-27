#n=0
#while read line;
#do
#    if [[ $line == "Cronbach"* ]]; then continue; fi
#    compare=$(echo $line'<'REPLACECRONBACHCUTOFF | bc)
#    if [[ compare -eq 1 ]]; then break; fi
#        ((n=$n+1))
#done < cronbach.txt

cat < eigenvectors.txt | cut -f1-100 > MetaBrain.eigenvectors.${n_PCs}PCs.txt
mv eigenvectors.txt $output_dir/7_evd_on_correlation_matrix/MetaBrain.eigenvectors.txt
mv pc-scores.txt $output_dir/7_evd_on_correlation_matrix/MetaBrain.pc-scores.txt
gzip $output_dir/7_evd_on_correlation_matrix/MetaBrain.eigenvectors.txt
gzip $output_dir/7_evd_on_correlation_matrix/MetaBrain.pc-scores.txt
# need unzippd file for next step, so only use zcat
zcat $output_dir/7_evd_on_correlation_matrix/MetaBrain.eigenvectors.${n_PCs}PCs.txt > $output_dir/7_PCA_on_correlation_matrix/MetaBrain.eigenvectors.${n_PCs}PCs.txt.gz
cd -
