
# This script runs 3rd step to making GeneNetwork: quantile normalisation before PCA

set -e
set -u

project_dir=
outdir=
jardir=
expression_file=
mem=
main(){
    module load Java/1.8.0_144-unlimited_JCE
    parse_commandline "$@"

    java -Xmx$mem -Xms$mem -jar $jardir/eqtl-mapping-pipeline.jar \
        --mode normalize \
        --in $expression_file \
        --out $outdir \
        --qqnorm \
        --logtransform \
        --adjustPCA \
        --maxnrpcaremoved 0 \
        --stepsizepcaremoval 0

    zcat $outdir/*.QuantileNormalized.Log2Transformed.PCAOverSamplesEigenvectors.txt.gz | awk 'BEGIN {OFS="\t"} NR == 1 {print "Sample","Comp1","Comp2"} NR > 1 {print $1,$2,$3}' > $outdir/pc1_2.txt

R --vanilla << EOF

    args = commandArgs(trailingOnly=TRUE)
    pcs <- read.table("$outdir/pc1_2.txt", header=T, sep="\t")
    png("$outdir/pca.png",width=300, height=300)
    plot(pcs\$Comp1, pcs\$Comp2)
    dev.off()

EOF
echo "PCA  wrtten to $outdir/pca.png"
}

usage(){
    # print the usage of the programme
    programname=$0
    echo "usage: $programname -e expression_file -p project_directory -o output_dir -m mem"
    echo "  -e      Expression file to remove duplciates from"
    echo "  -p      Base of the project_dir where config files will be written"
    echo "  -o      Output file that will be written"
    echo "  -j      Location of eqtl-mapping-pipeline.jar"
    echo "  -m      Memory to give to the eqtlgen jar file when running, will be appendended to -Xmx, so 8g will be java -jar -Xmx8g"
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
            -j | --jardir )             shift
                                        jardir=$1
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
    if [ -z "$jardir" ];
    then
        echo "ERROR: -j/--jardir not set!"
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






