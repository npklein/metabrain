
# This script runs 6th step to making GeneNetwork: Make the correlation matrix

set -e
set -u

project_dir=
outfile=
config_templates=
jardir=
eigenvector_file=
threads=
mem=
qos=
github_dir=
main(){
    module load Java/1.8.0_144-unlimited_JCE
    parse_commandline "$@"

    eigname=$(basename $eigenvector_file)
    eigname=${eigname%.gz}
    eigname=${eigname%.txt}
    rsync -vP $config_templates/6_CorrelationMatrix.json $project_dir/configs/9_CorrelateEigenvectors.$eigname.json
    sed -i "s;REPLACEEXPRFILE;$eigenvector_file;" $project_dir/configs/9_CorrelateEigenvectors.$eigname.json
    sed -i "s;REPLACEOUTFILE;$outfile;" $project_dir/configs/9_CorrelateEigenvectors.$eigname.json
    sed -i "s;REPLACETHREADS;$threads;" $project_dir/configs/9_CorrelateEigenvectors.$eigname.json
    mkdir -p $(dirname $outfile)

    echo "Starting sbatch with:"
    echo "#!/bin/bash
#SBATCH --job-name=correlation_$eigname
#SBATCH --output=$(dirname $outfile)/correlation_$eigname.out
#SBATCH --error=$(dirname $outfile)/correlation_$eigname.err
#SBATCH --time=05:59:59
#SBATCH --cpus-per-task $threads
#SBATCH --mem ${mem}b
#SBATCH --nodes 1
#SBATCH --qos $qos

ml Java;
java -Xmx$mem -Xms$mem -jar $github_dir/RunV13.jar $project_dir/configs/9_CorrelateEigenvectors.$eigname.json

if [ $? -eq 0 ];
then
    echo "success!"
    touch $(dirname $outfile)/correlation_$eigname.finished
else
    echo "error!"
    exit 1;
fi

" > $(dirname $outfile)/correlation_$eigname.sh

    echo "start sbatch with:"
    echo "sbatch $(dirname $outfile)/correlation_$eigname.sh"
    sbatch $(dirname $outfile)/correlation_$eigname.sh

    echo "sleep 2 minutes before checking if correlation is done"
    dt=$(date '+%d/%m/%Y %H:%M:%S');
    echo "$dt"
    sleep 120
    echo "done sleeping, check again"
    while [ ! -f $(dirname $outfile)/correlation_$eigname.finished ]
    do
      echo "$(dirname $outfile)/correlation_$eigname.finished does not exist yet"
      echo "sleep 2 minutes before checking again"
      dt=$(date '+%d/%m/%Y %H:%M:%S');
      echo "$dt"
      sleep 120
      echo "----------------------------------"
      echo "done sleeping, check again"
    done

    if [ ! -f $outfile ];
    then
        echo "ERROR: $outfile not made!"
        exit 1;
    fi
    echo "finished, gzip the output file (can take some time)"
    gzip $outfile
    echo "done gzipping"

}

usage(){
    # print the usage of the programme
    programname=$0
    echo "usage: $programname -e eigenvector_file -p project_directory -o output_dir -m mem -q qos"
    echo "  -e      Expression file to remove duplciates from"
    echo "  -p      Base of the project_dir where config files will be written"
    echo "  -o      Output file that will be written"
    echo "  -c      Dir with configuration template files"
    echo "  -j      Location of V13 jar file"
    echo "  -t      Number of threads"
    echo "  -m      mem to use"
    echo "  -q      qos to use"
    echo "  -g      github_dir"
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
            -e | --eigenvector_file )    shift
                                        eigenvector_file=$1
                                        ;;
            -o | --outfile )            shift
                                        outfile=$1
                                        ;;
            -c | --config_templates )   shift
                                        config_templates=$1
                                        ;;
            -t | --threads )            shift
                                        threads=$1
                                        ;;
            -j | --jardir )             shift
                                        jardir=$1
                                        ;;
            -q | --qos )                shift
                                        qos=$1
                                        ;;
            -m | --mem )                shift
                                        mem=$1
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
    if [ -z "$eigenvector_file" ];
    then
        echo "ERROR: -e/--eigenvector_file not set!"
        usage
        exit 1;
    fi
    if [ -z "$outfile" ];
    then
        echo "ERROR: -o/--outfile not set!"
        usage
        exit 1;
    fi
    if [ -z "$jardir" ];
    then
        echo "ERROR: -j/--jardir not set!"
        usage
        exit 1;
    fi
    if [ -z "$threads" ];
    then
        echo "ERROR: -t/--threads not set!"
        usage
        exit 1;
    fi
    if [ -z "$config_templates" ];
    then
        echo "ERROR: -c/--config_templates not set!"
        usage
        exit 1;
    fi
    if [ -z "$mem" ];
    then
        echo "ERROR: -m/--mem not set!"
        usage
        exit 1;
    fi
    if [ -z "$qos" ];
    then
        echo "ERROR: -q/--qos not set!"
        usage
        exit 1;
    fi
    if [ -z "$github_dir" ];
    then
        echo "ERROR: -g/--github_dir not set!"
        usage
        exit 1;
    fi
}

# [[ ${BASH_SOURCE[0]} = "$0" ]] -> main does not run when this script is sourced
# main "$@" -> Send the arguments to the main function (this way project flow can be at top)
# exit -> safeguard against the file being modified while it is interpreted
[[ ${BASH_SOURCE[0]} = "$0" ]] && main "$@"; exit;






