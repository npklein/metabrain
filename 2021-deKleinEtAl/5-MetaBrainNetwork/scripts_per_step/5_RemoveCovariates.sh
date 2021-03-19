
# This script runs 5th step to making GeneNetwork: Remove covariates from normalized data

set -e
set -u

project_dir=
outdir=
config_templates=
expression_file=
covar_matrix=
threads=
github_dir=
mem=
main(){
    module load Java/1.8.0_144-unlimited_JCE
    module load Python
    parse_commandline "$@"

    rsync -vP $config_templates/5_RemoveCovariates.json $project_dir/configs/
    head -n1 $expression_file | sed -e 's;\t;\n;g' > $outdir/samples.txt
    subset_covar="$outdir/$(basename ${covar_matrix%.*}).subset.txt"
    python $github_dir/table_scripts/subset_covar_matrix_remove_noVariance_rows.py $covar_matrix $outdir/samples.txt $subset_covar
    sed -i "s;REPLACEEXPRFILE;$expression_file;" $project_dir/configs/5_RemoveCovariates.json
    sed -i "s;REPLACEOUTDIR;$outdir/;" $project_dir/configs/5_RemoveCovariates.json
    sed -i "s;REPLCACECOVARMATRIX;$subset_covar;" $project_dir/configs/5_RemoveCovariates.json

    mkdir -p $(dirname $outdir)

    java -Xmx$mem -Xms$mem -jar $github_dir//RunV13.jar $project_dir/configs/5_RemoveCovariates.json


}

usage(){
    # print the usage of the programme
    programname=$0
    echo "usage: $programname -e expression_file -p project_directory -o output_dir -c config_dir -g github_dir -z covariance matrix -m memory"
    echo "  -e      Expression file to remove duplciates from"
    echo "  -p      Base of the project_dir where config files will be written"
    echo "  -o      Output directory that will be written"
    echo "  -c      Dir with configuration template files"
    echo "  -g      Github dir where RunV13.jar is located"
    echo "  -z      Covariance matrix to regress out"
    echo "  -m      Memory to reserve"
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
            -p | --project_dir )        shift
                                        project_dir=$1
                                        ;;
            -e | --expression_file )    shift
                                        expression_file=$1
                                        ;;
            -o | --outdir )            shift
                                        outdir=$1
                                        ;;
            -c | --config_templates )   shift
                                        config_templates=$1
                                        ;;
            -g | --github_dir )         shift
                                        github_dir=$1
                                        ;;
            -z | --covar_matrix )       shift
                                        covar_matrix=$1
                                        ;;
            -m | --mem )                shift
                                        mem=$1
                                        ;;
            -h | --help )               usage
                                        exit
                                        ;;
            * )                         echo "ERROR: Undexpected argument: $1"
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
    if [ -z "$outdir" ];
    then
        echo "ERROR: -o/--outdir not set!"
        usage
        exit 1;
    fi
    if [ -z "$github_dir" ];
    then
        echo "ERROR: -g/--github_dir not set!"
        usage
        exit 1;
    fi
    if [ -z "$config_templates" ];
    then
        echo "ERROR: -c/--config_templates not set!"
        usage
        exit 1;
    fi
    if [ -z "$covar_matrix" ];
    then
        echo "ERROR: -z/--covar_matrix not set!"
        usage
        exit 1;
    fi
    if [ -z "$mem" ];
    then
        echo "ERROR: -m/--mem not set!"
        usage
        exit 1;
    fi
}

# [[ ${BASH_SOURCE[0]} = "$0" ]] -> main does not run when this script is sourced
# main "$@" -> Send the arguments to the main function (this way project flow can be at top)
# exit -> safeguard against the file being modified while it is interpreted
[[ ${BASH_SOURCE[0]} = "$0" ]] && main "$@"; exit;






