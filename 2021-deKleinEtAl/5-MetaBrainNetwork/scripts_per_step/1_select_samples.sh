
# This script runs first step to making GeneNetwork: select samples

set -e
set -u

project_dir=
outdir=
config_templates=
sample_file=
expression_file=
github_dir=
main(){
    module load Java/1.8.0_144-unlimited_JCE
    parse_commandline "$@"

    echo "rsync -vP $config_templates/1_config_SelectSamples.json $project_dir/configs/"
    rsync -vP $config_templates/1_config_SelectSamples.json $project_dir/configs/

    sed -i "s;REPLACEEXPRFILE;$expression_file;" $project_dir/configs/1_config_SelectSamples.json
    sed -i "s;REPLACEOUTPUTDIR;$outdir;" $project_dir/configs/1_config_SelectSamples.json
    sed -i "s;REPLACEINCLUDESAMPLES;$sample_file;" $project_dir/configs/1_config_SelectSamples.json

    mkdir -p $(dirname $outdir)
    java -jar $github_dir/RunV13.jar $project_dir/configs/1_config_SelectSamples.json
}

usage(){
    # print the usage of the programme
    programname=$0
    echo "usage: $programname -e expression_file -p project_directory -o output_dir -s sample_file -g github_dir -c conf_templates"
    echo "  -e      Expression file to remove duplciates from"
    echo "  -p      Base of the project_dir where config files will be written"
    echo "  -o      Output file that will be written"
    echo "  -c      Dir with configuration template files"
    echo "  -g      location of the github directory brain_eQTL/"
    echo "  -s      File with samples to analyze"
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
            -g | --github_dir )             shift
                                        github_dir=$1
                                        ;;
            -s | --sample_file )        shift
                                        sample_file=$1
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
    if [ -z "$sample_file" ];
    then
        echo "ERROR: -s/--sample_file not set!"
        usage
        exit 1;
    fi
}

# [[ ${BASH_SOURCE[0]} = "$0" ]] -> main does not run when this script is sourced
# main "$@" -> Send the arguments to the main function (this way project flow can be at top)
# exit -> safeguard against the file being modified while it is interpreted
[[ ${BASH_SOURCE[0]} = "$0" ]] && main "$@"; exit;






