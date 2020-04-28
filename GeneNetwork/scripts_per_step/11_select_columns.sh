
set -p
set -u

project_dir=
outdir=
config_templates=
github_dir=
prediction_file=
auc_file=
main(){
    module load Java/1.8.0_144-unlimited_JCE
    parse_commandline "$@"
    rsync -vP $config_templates/11_config_MatrixScripts.GetCols_server.json $project_dir/configs/

    sed -i "s;REPLACEIN1;$prediction_file;" $project_dir/configs/11_config_MatrixScripts.GetCols_server.json
    sed -i "s;REPLACEIN2;$auc_file;" $project_dir/configs/11_config_MatrixScripts.GetCols_server.json
    sed -i "s;REPLACEOUT;$outdir/;" $project_dir/configs/11_config_MatrixScripts.GetCols_server.json

    mkdir -p $(dirname $outdir)

    java -jar $github_dir/RunV13.jar $project_dir/configs/11_config_MatrixScripts.GetCols_server.json


}

usage(){
    # print the usage of the programme
    programname=$0
    echo "usage: $programname -q prediction_file -p project_directory -o output_dir -c config_dir -a auc_file"
    echo "  -q      Prediction file to select signif only from"
    echo "  -p      Base of the project_dir where config files will be written"
    echo "  -o      Output directory that will be written"
    echo "  -c      Dir with configuration template files"
    echo "  -g      Location of V13 jar file"
    echo "  -a      AUC file"
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
            -q | --prediction_file )    shift
                                        prediction_file=$1
                                        ;;
            -a | --auc_file )           shift
                                        auc_file=$1
                                        ;;
            -o | --outdir )             shift
                                        outdir=$1
                                        ;;
            -c | --config_templates )   shift
                                        config_templates=$1
                                        ;;
            -g | --github_dir )         shift
                                        github_dir=$1
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
    if [ -z "$prediction_file" ];
    then
        echo "ERROR: -q/--prediction_file not set!"
        usage
        exit 1;
    fi
    if [ -z "$auc_file" ];
    then
        echo "ERROR: -a/--auc_file not set!"
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
}

# [[ ${BASH_SOURCE[0]} = "$0" ]] -> main does not run when this script is sourced
# main "$@" -> Send the arguments to the main function (this way project flow can be at top)
# exit -> safeguard against the file being modified while it is interpreted
echo $@
[[ ${BASH_SOURCE[0]} = "$0" ]] && main "$@"; exit;






