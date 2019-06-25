
# This runs the complete GeneNetwork processing pipeline

set -e
set -u

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
####

main(){
    parse_commandline "$@"
    print_command_arguments
    mkdir -p $output_dir
    1_select_samples
    2_remove_duplicate_samples
    3_quantileNormalized
    4_RemoveCovariates
}

print_command_arguments(){
    echo "Starting program with:"

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
}


1_select_samples(){
    # Step 1. Select samples from expression file
    output_dir_step1="${output_dir}/1_selectSamples/$(basename ${expression_file%.txt.gz}).extractedColumns.txt.gz"
    new_output_file_step1=$output_dir/1_selectSamples/$(basename ${expression_file%.txt.gz})_extractedColumnsnoVarianceRowsRemoved.txt.gz
    if [ ! -f ${new_output_file_step1} ]
    then
        echo "start step 1"
        bash $github_dir/GeneNetwork/scripts_per_step/1_select_samples.sh \
            -e $expression_file \
            -p $project_dir \
            -o $TMPDIR/1_selectSamples/ \
            -c $config_template_dir \
            -j $jar_dir \
            -s $sample_file
        # input of next file expect postfix of extractedColumnsnoVarianceRowsRemoved.txt.gz but is extractedColumns_noVarRemoved.txt.gz
        # change this
        f1=$TMPDIR/1_selectSamples/$(basename ${expression_file%.txt.gz})_extractedColumns_noVarRemoved.txt.gz
        mv $TMPDIR/1_selectSamples/ $output_dir/
        mv $f1 $TMPDIR/1_selectSamples/$(basename $new_output_file_step1)
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
            sbatch $TMPDIR/3_quantileNormalized/run_PCA.sh
            mv $MPDIR/3_quantileNormalized/ $output_dir/
        fi
    fi
}

4_RemoveCovariates(){
    output_file_step5=$output_dir/5_covariatesRemoved/$(basename ${output_file_step4%.txt.gz}.ProbesWithZeroVarianceRemoved.CovariatesRemoved.txt.gz)
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

            # The export.sh file has hardcoded paths for running PCA, change these
            rsync -vP $github_dir/GeneNetwork/scripts_per_step/run_PCA.sh $TMPDIR/3_quantileNormalized/run_PCA.sh
            sed -i "s;REPLACEGENECOVARIANCE;$TMPDIR/3_quantileNormalized/gene_covariance.txt;" $TMPDIR/3_quantileNormalized/run_PCA.sh
            sed -i "s;REPLACEOUT;$TMPDIR/3_quantileNormalized//;" $TMPDIR/3_quantileNormalized/run_PCA.sh
            sed -i "s;REPLACEPRECOR;$TMPDIR/3_quantileNormalized/pre_Correlation_Or_Covariance.txt;" $TMPDIR/3_quantileNormalized/run_PCA.sh
            sbatch $TMPDIR/3_quantileNormalized/run_PCA.sh
            mv $MPDIR/3_quantileNormalized/ $output_dir/
}

usage(){
    # print the usage of the programme
    programname=$0
    echo "usage: $programname -t TMPDIR -e expression_file -o output_dir -p project_dir -c config_template_dir -j jar_dir -s sample_file -g github_dir -z covar_table"
    echo "  -t      TMPDIR where files will be written during runtime"
    echo "  -e      Expression file"
    echo "  -p      Base of the project_dir where config files will be written"
    echo "  -c      Dir with configuration template files"
    echo "  -o      Output directory where results will be written"
    echo "  -j      Location of V13 jar file"
    echo "  -s      File with samples to include"
    echo "  -g      Github GeneNetwork directory"
    echo "  -z      Covariate table"
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
            -z | --covar_table )            shift
                                            covar_table=$1
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
}

# [[ ${BASH_SOURCE[0]} = "$0" ]] -> main does not run when this script is sourced
# main "$@" -> Send the arguments to the main function (this way project flow can be at top)
# exit -> safeguard against the file being modified while it is interpreted
[[ ${BASH_SOURCE[0]} = "$0" ]] && main "$@"; exit;






