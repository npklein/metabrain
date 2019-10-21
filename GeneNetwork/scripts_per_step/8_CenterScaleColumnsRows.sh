
# This script runs 8th step to making GeneNetwork: Center and scale first column wise then row wise the PC scores

set -e
set -u

infile=
outfile=
main(){
    module load RPlus/3.5.1-foss-2015b-v19.04.1
    parse_commandline "$@"

    file_dir=$(dirname "$0")
    mkdir -p $(dirname $outfile)
    Rscript $file_dir/../table_scripts/scale_columns_and_rows.R -i $infile -o ${outfile%.gz}
    gzip ${outfile%.gz}

    if [ $? -eq 0 ];
    then
        echo "success!"
    else
        echo "error!"
        exit 1;
    fi

}

usage(){
    # print the usage of the programme
    programname=$0
    echo "usage: $programname -i infile -o outfile"
    echo "  -i      PC scores, with genes on rows, PC scores on columns"
    echo "  -o      Output file that will be written"
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
            -i | --infile )        shift
                                        infile=$1
                                        ;;
            -o | --eigenvector_file )   shift
                                        outfile=$1
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
    if [ -z "$infile" ];
    then
        echo "ERROR: -i/--infile not set!"
        usage
        exit 1;
    fi
    if [ -z "$outfile" ];
    then
        echo "ERROR: -o/--outfile not set!"
        usage
        exit 1;
    fi
}

# [[ ${BASH_SOURCE[0]} = "$0" ]] -> main does not run when this script is sourced
# main "$@" -> Send the arguments to the main function (this way project flow can be at top)
# exit -> safeguard against the file being modified while it is interpreted
[[ ${BASH_SOURCE[0]} = "$0" ]] && main "$@"; exit;






