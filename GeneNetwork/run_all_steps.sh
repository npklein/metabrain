
# This runs the complete GeneNetwork processing pipeline

set -e
set -u

module load Python

#### Only lines that should be changed ####
# TMPDIR for writing files
TMPDIR=
# original expression table
expression_file=
# Directory to write output to
output_dir=
# Directory to write configuration files to (will be put in $project_dir/configs/)
project_dir=
# directory with config file templates
config_template_dir=
# directory containing V13.jar
jar_dir=
# file with samples to analyze
sample_file=
# github dir
github_dir=
# quant normalize only necesarry to select samples. we already know which are the good ones, so skip
quant_norm=false
# covariate table
covar_table=
# GTF file
gtf=
# number of threads to use for correlation and PCA step
threads=
# Cut-off value for cronbach alphas of PCA
cronbach_cutoff=
####

main(){
    parse_commandline "$@"
    print_command_arguments
    mkdir -p $output_dir
    echo "0_remove_genes"
    0_remove_genes
    echo "1_select_samples"
    1_select_samples
    echo "2_remove_duplicate_samples"
    2_remove_duplicate_samples
    echo "3_quantileNormalized"
    3_quantileNormalized
    echo "4_DeseqNormalizedData"
    4_DeseqNormalizedData
    echo "5_RemoveCovariates"
    5_RemoveCovariates
    echo "6_CorrelationMatrix"
    6_CorrelationMatrix
    echo "7_PCA_on_correlation"
    7_PCA_on_correlation
}


print_command_arguments(){
    echo "Processing expression data for GeneNetwork."
    echo "Starting program with:"
    echo "expression_file=$expression_file"
    echo "output_dir=$output_dir"
    echo "project_dir=$project_dir"
    echo "config_template_dir=$config_template_dir"
    echo "jar_dir=$jar_dir"
    echo "sample_file=$sample_file"
    echo "github_dir=$github_dir"
    echo "covar_table=$covar_table"
    echo "quant_norm=$quant_norm"
    echo "gtf=$gtf"
    echo "cronbach_cutoff=$cronbach_cutoff"
}


0_remove_genes(){
    # Step 1. Select samples from expression file
    outfile_step0="${output_dir}/0_remove_genes/$(basename ${expression_file%.txt.gz}).genesRemoved.txt.gz"
    if [ ! -f ${outfile_step0} ]
    then
        echo "start step 0"
        python $github_dir/GeneNetwork/scripts_per_step/0_remove_duplicate_sameCounts_scaffolded_genes.py \
            -e $expression_file \
            -o $TMPDIR/0_remove_genes/$(basename $outfile_step0) \
            -g $gtf
        mv $TMPDIR/0_remove_genes ${output_dir}/
    fi
}

1_select_samples(){
    # Step 1. Select samples from expression file
    new_output_file_step1=$output_dir/1_selectSamples/$(basename ${outfile_step0%.txt.gz})_extractedColumnsnoVarianceRowsRemoved.txt.gz
    echo $new_output_file_step1
    if [ ! -f ${new_output_file_step1} ]
    then
        echo "start step 1"
        bash $github_dir/GeneNetwork/scripts_per_step/1_select_samples.sh \
            -e $outfile_step0 \
            -p $project_dir \
            -o $TMPDIR/1_selectSamples/ \
            -c $config_template_dir \
            -j $jar_dir \
            -s $sample_file
        # input of next file expect postfix of extractedColumnsnoVarianceRowsRemoved.txt.gz but is extractedColumns_noVarRemoved.txt.gz
        # change this
        f1=$output_dir/1_selectSamples/$(basename ${outfile_step0%.txt.gz})_extractedColumns_noVarRemoved.txt.gz
        echo  "mv $TMPDIR/1_selectSamples/ $output_dir/"
        mv $TMPDIR/1_selectSamples/ $output_dir/
        mv $f1 $new_output_file_step1
        if [ ! -f ${new_output_file_step1} ];
        then
            echo "${new_output_file_step1} does not exist"
            exit 1;
        fi
    fi
}

2_remove_duplicate_samples(){
    # We have some ENA samples in our data, so run for duplicates. Add the relevant input files to the config file
    # Step 2. remove duplicates
    output_file_step2="${output_dir}/2_removeDuplicates/$(basename ${expression_file%.txt.gz}).duplicateSamplesRemoved_extractedColumns.txt.gz"
    if [ ! -f ${output_file_step2} ];
    then
        echo "start step 2"
        bash $github_dir/GeneNetwork/scripts_per_step/2_remove_duplicate_samples.sh \
            -p $project_dir \
            -e $new_output_file_step1 \
            -o $TMPDIR/$(basename $output_file_step2) \
            -c $config_template_dir \
            -j $jar_dir
        echo "done"
        mkdir -p $(dirname $output_file_step2)
        mv $TMPDIR/$(basename $output_file_step2) $output_file_step2
    fi
}

3_quantileNormalized(){
    if [ "quant_norm" = true ];
    then
        # Step 3. make PCA of quantile normalized data. Do visual inspection, if samples need to be removed, run from step 2 again (selecting only correct samples)
        # and then skip step 1 and 3 and go to step 4. I think this step is nececarry for step 4, not sure why though.
        output_file_step3=$output_dir/3_quantileNormalized/$(basename ${output_file_step2%.txt.gz}.QuantileNormalized.txt.gz)
        if [ ! -f $output_file_step3 ];
        then
            bash $github_dir/GeneNetwork/scripts_per_step/3_PCA_on_quantNormalizedData.sh \
                -p $project_dir \
                -e $output_file_step2 \
                -o $TMPDIR/3_quantileNormalized/ \
                -c $config_template_dir \
                -j $jar_dir

            # The export.sh file has hardcoded paths for running PCA, change these
            rsync -vP $github_dir/GeneNetwork/scripts_per_step/run_PCA.sh $TMPDIR/3_quantileNormalized/run_PCA.sh
            sed -i "s;REPLACEGENECOVARIANCE;$TMPDIR/3_quantileNormalized/gene_covariance.txt;" $TMPDIR/3_quantileNormalized/run_PCA.sh
            sed -i "s;REPLACEOUT;$TMPDIR/3_quantileNormalized//;" $TMPDIR/3_quantileNormalized/run_PCA.sh
            sed -i "s;REPLACEPRECOR;$TMPDIR/3_quantileNormalized/pre_Correlation_Or_Covariance.txt;" $TMPDIR/3_quantileNormalized/run_PCA.sh
            sbatch --wait $TMPDIR/3_quantileNormalized/run_PCA.sh
            wait
            mv $MPDIR/3_quantileNormalized/ $output_dir/
        fi
    fi
}

4_DeseqNormalizedData(){
    output_file_step4="$output_dir/4_deseqNormalized/$(basename ${output_file_step2%.txt.gz}noVarianceRowsRemoved.DESeqNorm.log2.txt)"
    if [ ! -f $output_file_step4 ];
    then
        # Step 4. Do deseq normalisation on the original counts
        bash $github_dir/GeneNetwork/scripts_per_step/4_DeseqNormalizedData.sh \
            -p $project_dir \
            -e $output_file_step2 \
            -f $TMPDIR/4_deseqNormalized/$(basename ${output_file_step2%.txt.gz}noVarianceRowsRemoved.DESeqNorm.txt.gz) \
            -o $TMPDIR/4_deseqNormalized/ \
            -c $config_template_dir \
            -j $jar_dir \
            -z $TMPDIR/4_deseqNormalized/PCA_corrected_expression/

        mv $TMPDIR/4_deseqNormalized/ $output_dir/
    fi
}

5_RemoveCovariates(){
    output_file_step5=$output_dir/5_covariatesRemoved/$(basename ${output_file_step4%.txt}.ProbesWithZeroVarianceRemoved.CovariatesRemoved.txt.gz)
    echo "$output_file_step5"
    if [ ! -f $output_file_step5 ];
    then
        # Step 5. Remove covariates from deseq normalized data
        bash $github_dir/GeneNetwork/scripts_per_step/5_RemoveCovariates.sh \
            -p $project_dir \
            -e $output_file_step4 \
            -o $TMPDIR/5_covariatesRemoved \
            -c $config_template_dir \
            -j $jar_dir \
            -z $covar_table
        mv $TMPDIR/5_covariatesRemoved $output_dir/
    fi
}

6_CorrelationMatrix(){
    output_file_step6="$output_dir/6_correlation_matrix/MetaBrain.deseqNorm.covarCorrected.correlation.txt"
    if [ ! -f $output_file_step6 ]  && [ ! -f ${output_file_step6}.gz ];
    then
        # Step 6. Make correlation matrix
        bash $github_dir/GeneNetwork/scripts_per_step/6_CorrelationMatrix.sh \
            -p $project_dir \
            -e $output_file_step5 \
            -o ${output_file_step6} \
            -c $config_template_dir \
            -j $jar_dir \
            -t $threads
    fi
}

7_PCA_on_correlation(){
    if [ ! -f $output_file_step6 ] && [ -f ${output_file_step6}.gz ];
    then
        zcat ${output_file_step6}.gz > ${output_file_step6}
    fi

    if [ ! -f ${output_file_step5%.gz} ] && [ -f ${output_file_step5} ];
    then
         zcat ${output_file_step5} > ${output_file_step5%.gz}
    fi


    if [ ! -f $output_dir/7_PCA_on_correlation_matrix/MetaBrain.pc-scores.txt ];
    then
        # step 7. Run PCA on correlation matrix
        mkdir -p  $output_dir/7_PCA_on_correlation_matrix/
        rsync -vP $github_dir/GeneNetwork/scripts_per_step/7_PCA_on_correlation.sh $output_dir/7_PCA_on_correlation_matrix/7_PCA_on_correlation.sh
        echo "replace"
        sed -i "s;REPLACEPCASETUP;$github_dir/GeneNetwork/PCA_setup_calculon.sh;" $output_dir/7_PCA_on_correlation_matrix/7_PCA_on_correlation.sh
        sed -i "s;REPLACECORMATRIX;$output_file_step6;" $output_dir/7_PCA_on_correlation_matrix/7_PCA_on_correlation.sh
        sed -i "s;REPLACECOVCORRECTEDEXPRESSION;${output_file_step5%.gz};" $output_dir/7_PCA_on_correlation_matrix/7_PCA_on_correlation.sh
        sed -i "s;REPLACEOUTDIR;$output_dir/;" $output_dir/7_PCA_on_correlation_matrix/7_PCA_on_correlation.sh
        sed -i "s;REPLACETHREADS;$threads;" $output_dir/7_PCA_on_correlation_matrix/7_PCA_on_correlation.sh
        sed -i "s;REPLACECRONBACHCUTOFF;${cronbach_cutoff};" $output_dir/7_PCA_on_correlation_matrix/7_PCA_on_correlation.sh
        echo "!!!PCA++ does not work correctly on the cluster. Copy data to GeneNetwork and run 7_PCA_on_correlation.sh!!!"
        echo "Afterwards, copy data back and run run_all_steps.sh again"
        exit 0;
        cd $output_dir/7_PCA_on_correlation_matrix/
        sbatch 7_PCA_on_correlation.sh
        cd -
    fi
}


usage(){
    # print the usage of the programme
    programname=$0
    echo "usage: $programname -t TMPDIR -e expression_file -o output_dir -p project_dir -c config_template_dir -j jar_dir -s sample_file -g github_dir -z covar_table -t threads -a cronbach_cutoff -q quant_norm (default: false)"
    echo "  -t      TMPDIR where files will be written during runtime"
    echo "  -e      Expression file"
    echo "  -p      Base of the project_dir where config files will be written"
    echo "  -c      Dir with configuration template files"
    echo "  -o      Output directory where results will be written"
    echo "  -j      Location of V13 jar file"
    echo "  -s      File with samples to include"
    echo "  -g      Github GeneNetwork directory"
    echo "  -z      Covariate table"
    echo "  -a      GTF file"
    echo "  -v      Number of threads to use for correlation step and PCA step"
    echo "  -q      true/false wether quntile normalization should be done"
    echo "  -a      Cronbach alpha cut-off"
    echo "  -h      display help"
    exit 1
}

parse_commandline(){
    # Check to see if at least one argument is given
    if [ $# -eq 0 ]
    then
        echo "ERROR: No arguments supplied"
        usage
        exit 1;
    fi

    while [[ $# -ge 1 ]]; do
        case $1 in
            -t | --TMPDIR )                 shift
                                            TMPDIR=$1
                                            ;;
            -p | --project_dir )            shift
                                            project_dir=$1
                                            ;;
            -e | --expression_file )        shift
                                            expression_file=$1
                                            ;;
            -o | --output_dir )             shift
                                            output_dir=$1
                                            ;;
            -c | --config_template_dir )    shift
                                            config_template_dir=$1
                                            ;;
            -j | --jar_dir )                shift
                                            jar_dir=$1
                                            ;;
            -s | --sample_file )            shift
                                            sample_file=$1
                                            ;;
            -g | --github_dir )             shift
                                            github_dir=$1
                                            ;;
            -a | --gtf )                    shift
                                            gtf=$1
                                            ;;
            -v | --threads )                shift
                                            threads=$1
                                            ;;
            -z | --covar_table )            shift
                                            covar_table=$1
                                            ;;
            -a | --cronbach_cutoff )        shift
                                            cronbach_cutoff=$1
                                            ;;
            -q | --quant_norm )             shift
                                            quant_norm=$1
                                            ;;
            -h | --help )                   usage
                                            exit
                                            ;;
            * )                             echo "ERROR: Undexpected argument: $1"
                                            usage
                                            exit 1
        esac
        shift
    done

    # if -z tests if variable is empty. Make sure the relevant variables are set
    if [ -z "$project_dir" ];
    then
        echo "ERROR: -p/--project_dir not set!"
        usage
        exit 1;
    fi
    if [ -z "$expression_file" ];
    then
        echo "ERROR: -e/--expression_file not set!"
        usage
        exit 1;
    fi
    if [ -z "$output_dir" ];
    then
        echo "ERROR: -o/--output_dir not set!"
        usage
        exit 1;
    fi
    if [ -z "$jar_dir" ];
    then
        echo "ERROR: -j/--jar_dir not set!"
        usage
        exit 1;
    fi
    if [ -z "$config_template_dir" ];
    then
        echo "ERROR: -c/--config_template_dir not set!"
        usage
        exit 1;
    fi
        if [ -z "$sample_file" ];
    then
        echo "ERROR: -s/--sample_file not set!"
        usage
        exit 1;
    fi
    if [ -z "$github_dir" ];
    then
        echo "ERROR: -g/--github_dir not set!"
        usage
        exit 1;
    fi
    if [ -z "$gtf" ];
    then
        echo "ERROR: -a/--gtf not set!"
        usage
        exit 1;
    fi
    if [ -z "$config_template_dir" ];
    then
        echo "ERROR: -c/--config_template_dir not set!"
        usage
        exit 1;
    fi
    if [ -z "$covar_table" ];
    then
        echo "ERROR: -z/--covar_table not set!"
        usage
        exit 1;
    fi
    if [ -z "$crobach_cutoff" ];
    then
        echo "ERROR: -a/--cronbach_cutoff not set!"
        usage
        exit 1;
    fi
    if [ -z "$threads" ];
    then
        echo "ERROR: -v/--threads not set!"
        usage
        exit 1;
    fi
}

# [[ ${BASH_SOURCE[0]} = "$0" ]] -> main does not run when this script is sourced
# main "$@" -> Send the arguments to the main function (this way project flow can be at top)
# exit -> safeguard against the file being modified while it is interpreted
[[ ${BASH_SOURCE[0]} = "$0" ]] && main "$@"; exit;






