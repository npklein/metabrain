set -e
set -u
#### Only lines that should be changed ####
# TMPDIR for writing files
TMPDIR=/groups/umcg-pub/scr01/
# original expression table
input_expression_file="/groups/umcg-pub/tmp04/GeneNetworkTest/data/ten_ENA_samples.txt.gz"
# Directory to write output to
output_dir="/groups/umcg-pub/tmp04/GeneNetworkTest/output/"
# Directory to write configuration files to (will be put in $project_dir/configs/)
project_dir="/groups/umcg-pub/tmp04/GeneNetworkTest/"
# directory with config file templates
config_template_dir="/groups/umcg-pub/tmp04/GeneNetworkTest/brain_eQTL/GeneNetwork/config_file_templates/"
# directory containing V13.jar
jar_dir="/groups/umcg-pub/tmp04/GeneNetworkTest/jars"
# file with samples to analyze
sample_file="/groups/umcg-pub/tmp04/GeneNetworkTest/data/samples_to_include.txt"
# github dir
github_dir=/groups/umcg-pub/tmp04/GeneNetworkTest/brain_eQTL/
# quant normalize only necesarry to select samples. we already know which are the good ones, so skip
quant_norm=false
# covariate table
covar_table="/groups/umcg-biogen/tmp04/umcg-ndeklein/GeneNetwork/data/2019-04-11-Freeze2.TMM.Covariates-Numeric-Top10Covariates.transposed.txt"
####

mkdir -p $output_dir

# We have some ENA samples in our data, so run for duplicates. Add the relevant input files to the config file
# Step 1. remove duplicates
outfile_step1="${output_dir}/1_removeDuplicates/$(basename ${input_expression_file%.txt.gz}).duplicateSamplesRemoved.txt.gz"
if [ ! -f ${outfile_step1} ];
then
    echo "start step 1"
    bash $github_dir/GeneNetwork/scripts_per_step/1_remove_duplicate_samples.sh \
        -p $project_dir \
        -e $input_expression_file \
        -o $TMPDIR/$(basename $outfile_step1) \
        -c $config_template_dir \
        -j $jar_dir
    echo "done"
    echo ""
    echo ""
    echo "Do: mv $TMPDIR/$(basename $outfile_step1) $outfile_step1"
    mkdir -p $(dirname $outfile_step1)
    mv $TMPDIR/$(basename $outfile_step1) $outfile_step1
fi

# Step 2. Select samples from expression file
outfile_step2="${output_dir}/2_selectSamples/$(basename ${input_expression_file%.txt.gz}).duplicateSamplesRemoved_extractedColumns.txt.gz"
new_outfile_step2=$output_dir/2_selectSamples/$(basename ${input_expression_file%.txt.gz}).duplicateSamplesRemoved_extractedColumnsnoVarianceRowsRemoved.txt.gz
if [ ! -f ${new_outfile_step2} ]
then
    echo "start step 2"
    bash $github_dir/GeneNetwork/scripts_per_step/2_select_samples.sh \
        -p $project_dir \
        -e $outfile_step1 \
        -o $TMPDIR/2_selectSamples/ \
        -c $config_template_dir \
        -j $jar_dir \
        -s $sample_file
    # input of next file expect postfix of extractedColumnsnoVarianceRowsRemoved.txt.gz but is extractedColumns_noVarRemoved.txt.gz
    # change this
    f1=$TMPDIR/2_selectSamples/$(basename ${input_expression_file%.txt.gz}).duplicateSamplesRemoved_extractedColumns_noVarRemoved.txt.gz
    mv $f1 $TMPDIR/2_selectSamples/$(basename $new_outfile_step2)
    mv $TMPDIR/2_selectSamples/ $output_dir/
fi

if [ "quant_norm" = true ];
then
    # Step 3. make PCA of quantile normalized data. Do visual inspection, if samples need to be removed, run from step 2 again (selecting only correct samples)
    # and then skip step 1 and 3 and go to step 4. I think this step is nececarry for step 4, not sure why though.
    outfile_step3=$output_dir/3_quantileNormalized/$(basename ${new_outfile_step2%.txt.gz}.QuantileNormalized.txt.gz)
    if [ ! -f $outfile_step3 ];
    then
        bash $github_dir/GeneNetwork/scripts_per_step/3_PCA_on_quantNormalizedData.sh \
            -p $project_dir \
            -e $new_outfile_step2 \
            -o $TMPDIR/3_quantileNormalized/ \
            -c $config_template_dir \
            -j $jar_dir

        # The export.sh file has hardcoded paths for running PCA, change these
        rsync -vP $github_dir/GeneNetwork/scripts_per_step/run_PCA.sh $TMPDIR/3_quantileNormalized/3_run_PCA.sh
        sed -i "s;REPLACEGENECOVARIANCE;$TMPDIR/3_quantileNormalized/gene_covariance.txt;" $TMPDIR/3_quantileNormalized/3_run_PCA.sh
        sed -i "s;REPLACEOUT;$TMPDIR/3_quantileNormalized//;" $TMPDIR/3_quantileNormalized/3_run_PCA.sh
        sed -i "s;REPLACEPRECOR;$TMPDIR/3_quantileNormalized/pre_Correlation_Or_Covariance.txt;" $TMPDIR/3_quantileNormalized/3_run_PCA.sh
        bash $TMPDIR/3_quantileNormalized/3_run_PCA.sh
        mv $MPDIR/3_quantileNormalized/ $output_dir/
    fi
fi

outfile_step4=$output_dir/4_deseqNormalized/$(basename ${new_outfile_step2%.txt.gz}noVarianceRowsRemoved.DESeqNorm.txt.gz)
if [ ! -f $outfile_step4 ];
then
    # Step 4. Do deseq normalisation on the original counts
    bash $github_dir/GeneNetwork/scripts_per_step/4_DeseqNormalizedData.sh \
        -p $project_dir \
        -e $new_outfile_step2 \
        -o $TMPDIR/4_deseqNormalized/ \
        -c $config_template_dir \
        -j $jar_dir \
        -z $TMPDIR/4_deseqNormalized/PCA_corrected_expression/

    mv $TMPDIR/4_deseqNormalized/ $output_dir/
fi

outfile_step5=$output_dir/5_covariatesRemoved/$(basename ${outfile_step4%.txt.gz}.CovariatesRemoved.txt.gz)
if [ ! -f $outfile_step5 ];
then
    # Step 5. Remove covariates from deseq normalized data
    bash $github_dir/GeneNetwork/scripts_per_step/5_RemoveCovariates.sh \
        -p $project_dir \
        -e $outfile_step4 \
        -o $TMPDIR/5_covariatesRemoved \
        -c $config_template_dir \
        -j $jar_dir \
        -z $covar_table
    mv $TMPDIR/5_covariatesRemoved $output_dir/
fi

outfile_step6="$output_dir/6_correlation_matrix/MetaBrain.deseqNorm.covarCorrected.correlation.txt.gz"
if [ ! -f $outfile_step6 ];
then
    # Step 6. Make correlation matrix
    bash $github_dir/GeneNetwork/scripts_per_step/6_CorrelationMatrix.sh \
        -p $project_dir \
        -e $outfile_step5 \
        -o $TMPDIR/$(basename $outfile_step6) \
        -c $config_template_dir \
        -j $jar_dir
    mkdir -p $(dirname $outfile_step6)
    mv $TMPDIR/$(basename $outfile_step6) $outfile_step6
fi

if [ ! -f $output_dir/7_PCA_on_correlation_matrix/MetaBrain.pc-scores.txt ];
then
    # step 7. Run PCA on correlation matrix
    rsync -vP $github_dir/GeneNetwork/scripts_per_step/7_PCA_on_correlation.sh $TMPDIR/7_PCA_on_correlation.sh
    sed -i "s;REPLACEOUTDIR;$output_dir/;" $TMPDIR/7_PCA_on_correlation.sh
    sed -i "s;REPLACEOUTFILE;$outfile_step6;" $TMPDIR/7_PCA_on_correlation.sh
    sed -i "s;REPLACEPCASETUP;$github_dir/GeneNetwork/PCA_setup_calculon.sh;" $TMPDIR/7_PCA_on_correlation.sh
    sbatch $TMPDIR/7_PCA_on_correlation.sh
fi


