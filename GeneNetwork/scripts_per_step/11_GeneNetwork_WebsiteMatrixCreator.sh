
# This script runs 5th step to making GeneNetwork: Remove covariates from normalized data

set -e
set -u

project_dir=
outdir=
config_templates=
jardir=
identity_matrix=
term_file=
auc_file=
zscore_matrix=
main(){
    module load Java/1.8.0_144-unlimited_JCE
    parse_commandline "$@"

    rsync -vP $config_templates/10_config_GeneNetwork.WebsiteMatrixCreator_Server.json $project_dir/configs/

    sed -i "s;REPLACETERMFILE;$term_file;" $project_dir/configs/10_config_GeneNetwork.WebsiteMatrixCreator_Server.json
    sed -i "s;REPLACEIDENTITY;$identity_matrix;" $project_dir/configs/10_config_GeneNetwork.WebsiteMatrixCreator_Server.json
    sed -i "s;REPLACEAUC;$auc_file;" $project_dir/configs/10_config_GeneNetwork.WebsiteMatrixCreator_Server.json
    sed -i "s;REPLACEZSCORE;$zscore_matrix;" $project_dir/configs/10_config_GeneNetwork.WebsiteMatrixCreator_Server.json
    sed -i "s;REPLACEOUTDIR;$outdir;" $project_dir/configs/10_config_GeneNetwork.WebsiteMatrixCreator_Server.json

    mkdir -p $(dirname $outdir)

    java -jar $jardir/RunV13.jar $project_dir/configs/9_config_MatrixScripts.GetCols_server.json


}

usage(){
    # print the usage of the programme
    programname=$0
    echo "usage: $programname -i identity_matrix -t term_file -a auc_file -z zscore_matrix -p project_directory -o output_dir -j jardir"
    echo "  -i      Identitiy matrix"
    echo "  -p      Base of the project_dir where config files will be written"
    echo "  -o      Output directory that will be written"
    echo "  -c      Dir with configuration template files"
    echo "  -j      Location of V13 jar file"
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
            -i | --identity_matrix  )   shift
                                        identity_matrix=$1
                                        ;;
            -t | --term_file  )         shift
                                        term_file=$1
                                        ;;
            -a | --auc_file  )          shift
                                        auc_file=$1
                                        ;;
            -z | --zscore_matrix  )     shift
                                        zscore_matrix=$1
                                        ;;
            -o | --outdir )             shift
                                        outdir=$1
                                        ;;
            -c | --config_templates )   shift
                                        config_templates=$1
                                        ;;
            -j | --jardir )             shift
                                        jardir=$1
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
    if [ -z "$identity_matrix" ];
    then
        echo "ERROR: -i/--identity_matrix not set!"
        usage
        exit 1;
    fi
    if [ -z "$term_file" ];
    then
        echo "ERROR: -t/--term_file not set!"
        usage
        exit 1;
    fi
    if [ -z "$auc_file" ];
    then
        echo "ERROR: -a/--auc_file not set!"
        usage
        exit 1;
    fi
    if [ -z "$zscore_matrix" ];
    then
        echo "ERROR: -z/--zscore_matrix not set!"
        usage
        exit 1;
    fi
    if [ -z "$outdir" ];
    then
        echo "ERROR: -o/--outdir not set!"
        usage
        exit 1;
    fi
    if [ -z "$jardir" ];
    then
        echo "ERROR: -j/--jardir not set!"
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
[[ ${BASH_SOURCE[0]} = "$0" ]] && main "$@"; exit;






