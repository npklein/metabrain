#!/bin/bash
# This runs the complete GeneNetwork processing pipeline

set -e
set -u

module load Python

#### Only lines that should be changed ####
# TMPDIR for writing files
TMPDIR=
# Directory to write output to
output_dir=
# Directory to write configuration files to (will be put in $project_dir/configs/)
project_dir=
# directory containing V13.jar
jar_dir=
# github dir
github_dir=
# GTF file
gtf=
# number of threads to use for correlation and evd step
threads=
# gene_network_dir that contains GeneNetworkBackend-1.0.7-SNAPSHOT-jar-with-dependencies.jar and PathwayMatrix/
gene_network_dir=
# memory to use when doing normalization, evd etc
mem=
# Name of the project (e.g. GeneNetwork, MetaBrain, KidneyNetwork, etc)
# This needs to be same as previous step
name=
# Number of eigenvectors to keep
n_eigenvectors=
# qos when submitting to cluster
qos=regular
####

main(){
    parse_commandline "$@"
    print_command_arguments
    cd $github_dir/GeneNetwork
    if [ ! -f RunV13.jar ];
    then
        echo "RunV13.jar does not exist, download now"
        wget https://github.com/npklein/brain_eQTL/releases/download/2020-03-02/RunV13.jar
    fi
    md5=$(md5sum RunV13.jar | awk '{print $1}')
    if [ $md5 != "9b292956792206b7cdb8a4e7fee5a06b" ];
    then
        echo "ERROR: md5sum of RunV13.jar not correct"
        exit 1;
    fi
    cd -
    mkdir -p $output_dir
    echo "8_center_scale"
    8_center_scale
    echo "9_correlate_eigenvectors"
    9_correlate_eigenvectors
    echo "10_GeneNetwork_predictions"
    10_GeneNetwork_predictions
    echo "11_bonferonni_correction"
    11_bonferonni_correction
    echo "12_WebsiteMatrixCreator"
    12_WebsiteMatrixCreator
}


print_command_arguments(){
    echo "-------------------------------------------"
    echo "Processing expression data for GeneNetwork."
    echo "-------------------------------------------"
    echo "Starting program with:"
    echo "output_dir=$output_dir"
    echo "project_dir=$project_dir"
    echo "jar_dir=$jar_dir"
    echo "github_dir=$github_dir"
    echo "gtf=$gtf"
    echo "gene_network_dir=$gene_network_dir"
    echo "mem=$mem"
    echo "name=$name"
    echo "qos=$qos"
    echo "n_eigenvectors=$n_eigenvectors"
    echo "-------------------------------------------"
    echo ""
    echo ""
}

8_center_scale(){
    output_file_step8="$output_dir/8_CenterScaledColumnsRows/$name.${n_eigenvectors}_eigenvectors.columnsRowsCenterScaled.txt.gz"
    if [ ! -f ${output_file_step8} ];
    then
        echo "${output_file_step8}.gz does not exist yet, start step 8"
        if [ ! -f ${output_dir}/7_evd_on_correlation_matrix/$name.eigenvectors.${n_eigenvectors}_eigenvectors.txt.gz ];
        then
            nforcut=$(expr $n_eigenvectors + 1)
            zcat < $output_dir/7_evd_on_correlation_matrix/$name.eigenvectors.txt.gz | cut -f1-$nforcut > $output_dir/7_evd_on_correlation_matrix/$name.eigenvectors.${n_eigenvectors}_eigenvectors.txt
            gzip $output_dir/7_evd_on_correlation_matrix/$name.eigenvectors.${n_eigenvectors}_eigenvectors.txt
        fi
        bash $github_dir/GeneNetwork/scripts_per_step/8_CenterScaleColumnsRows.sh \
            -i $output_dir/7_evd_on_correlation_matrix/$name.eigenvectors.${n_eigenvectors}_eigenvectors.txt \
            -o ${output_file_step8}
    else
        echo "${output_file_step8}.gz exists, go to step 9"
    fi
}

9_correlate_eigenvectors(){
    output_file_step9="$output_dir/9_eigenvector_correlation_matrix/$name.eigenvectors.columnsRowsCenterScaled.correlation.txt"
    if [ ! -f ${output_file_step9}.gz ];
    then
        echo "${output_file_step9}.gz does not exist yet, start step 9"
        # Step 8. Make correlation matrix of eigenvectors
        bash $github_dir/GeneNetwork/scripts_per_step/9_CorrelateEigenvectors.sh \
            -p $project_dir \
            -e $output_file_step8 \
            -o ${output_file_step9} \
            -c $github_dir/GeneNetwork/config_file_templates/ \
            -j $jar_dir \
            -t $threads \
            -q $qos \
            -m $mem \
            -g $github_dir/GeneNetwork/
    else
        echo "${output_file_step9}.gz already exists, go to step 10"
    fi
}

10_GeneNetwork_predictions(){
    mkdir -p  $output_dir/10_GeneNetwork_predictions/scripts/
    declare -A fullname=( ["go_F"]="goa_human.gaf_F" ["go_P"]="goa_human.gaf_P" ["kegg"]="c2.cp.kegg.v6.1.entrez.gmt" ["hpo"]="ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt" ["reactome"]="Ensembl2Reactome_All_Levels.txt" ["go_C"]="goa_human.gaf_C")
    for type in "go_F" "go_P" "kegg" "hpo" "reactome" "go_C";
    do
        mkdir -p $output_dir/10_GeneNetwork_predictions/scripts/$type/
        outfile="$output_dir/10_GeneNetwork_predictions/$name.${type}_predictions.txt"
        if [ ! -f $outfile ];
        then
            for f in ${gene_network_dir}/PathwayMatrix/split_matrices/$type/*txt;
            do
                outfile="$output_dir/10_GeneNetwork_predictions/${type}/$(basename ${f%.txt}).predictions.txt"
                rsync -vP $github_dir/GeneNetwork/scripts_per_step/10_GeneNetwork_predictions.sh $output_dir/10_GeneNetwork_predictions/scripts/$type/${type}_$(basename ${f%.txt}).predictions.sh
                sed -i "s;REPLACENAME;${type}_$(basename ${f%.txt});" $output_dir/10_GeneNetwork_predictions/scripts/$type/${type}_$(basename ${f%.txt}).predictions.sh
                sed -i "s;REPLACEGENENETWORKDIR;${gene_network_dir};" $output_dir/10_GeneNetwork_predictions/scripts/$type/${type}_$(basename ${f%.txt}).predictions.sh
                sed -i "s;REPLACEIDENTITYMATRIX;$f;" $output_dir/10_GeneNetwork_predictions/scripts/$type/${type}_$(basename ${f%.txt}).predictions.sh
                sed -i "s;REPLACEEIGENVECTORS;${output_dir}/7_PCA_on_correlation_matrix/$name.eigenvectors.filteredOnPathwayGenes.txt;" $output_dir/10_GeneNetwork_predictions/scripts/$type/${type}_$(basename ${f%.txt}).predictions.sh
                sed -i "s;REPLACEOUT;$outfile;" $output_dir/10_GeneNetwork_predictions/scripts/$type/${type}_$(basename ${f%.txt}).predictions.sh
                sed -i "s;REPLACETYPE;${fullname[$type]};" $output_dir/10_GeneNetwork_predictions/scripts/$type/${type}_$(basename ${f%.txt}).predictions.sh
                cd $output_dir/10_GeneNetwork_predictions/scripts/$type/
                jobid=$(sbatch --parsable ${type}_$(basename ${f%.txt}).predictions.sh)
                cd -
            done
            for f in ${gene_network_dir}/PathwayMatrix/split_matrices/$type/*txt;
            do
                startDate=`date`
                echo "Start testing if ${type}_$(basename ${f%.txt}).sh is finished"
                cd $output_dir/10_GeneNetwork_predictions/scripts/$type/
                while [ ! $(squeue -j $jobid | wc -l) -eq 0 ];
                do
                    echo "sbatch starting time: $startDate"
                    echo "current time: `date`"
                    echo "Job still running. Sleep 2 minutes"
                    sleep 120
                done
                echo "Finished!"
                if [ ! -f ${type}_$(basename ${f%.txt}).finished ];
                then
                    echo "ERROR: check job error log what went wrong"
                    exit 1;
                fi
                grep AUC ${type}_$(basename ${f%.txt}).out | cut -f2,7,11,13 >> ${outfile%.txt}.AUC.txt;
                cd -
            done
        fi
    done
}

11_bonferonni_correction(){
    echo "!!!not implemented correctly yet!!!"
    exit 1;
}

11_WebsiteMatrixCreator(){
    echo "!!!not implemented correctly yet!!!"
    exit 1;
}

usage(){
    # print the usage of the programme
    programname=$0
    echo "usage: $programname -t TMPDIR -o output_dir -p project_dir"
    echo "                    -j jar_dir -g github_dir -v threads"
    echo "                    -d GeneNetworkDir -n name -m mem -z n_eigenvectors [-q qos]"
    echo "  -t      TMPDIR where files will be written during runtime"
    echo "  -e      Expression file"
    echo "  -p      Base of the project_dir where config files will be written"
    echo "  -o      Output directory where results will be written"
    echo "  -j      Location of eqtl-mapping-pipeline.jar"
    echo "  -s      File with samples to include"
    echo "  -g      Location of git cloned https://github.com/npklein/brain_eQTL/ directory"
    echo "  -a      GTF file"
    echo "  -v      Number of threads to use for correlation step and PCA step"
    echo "  -d      GeneNetwork directory (with backend data for predictions"
    echo "  -m      Memory to use for some steps"
    echo "  -n      Name that will be used in output file. Needs to be same as for the previous script!"
    echo "  -z      Number of eigenvectors to use"
    echo "  -q      qos to use when submitting to cluster"
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
            -o | --output_dir )             shift
                                            output_dir=$1
                                            ;;
            -j | --jar_dir )                shift
                                            jar_dir=$1
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
            -d | --gene_network_dir )       shift
                                            gene_network_dir=$1
                                            ;;
            -m | --mem )                    shift
                                            mem=$1
                                            ;;
            -n | --name )                   shift
                                            name=$1
                                            ;;
            -z | --n_eigenvectors )         shift
                                            n_eigenvectors=$1
                                            ;;
            -q | --qos )                    shift
                                            qos=$1
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
    if [ -z "$gene_network_dir" ];
    then
        echo "ERROR: -d/--gene_network_dir not set!"
        usage
        exit 1;
    fi
    if [ -z "$threads" ];
    then
        echo "ERROR: -v/--threads not set!"
        usage
        exit 1;
    fi
    if [ -z "$mem" ];
    then
        echo "ERROR: -m/--mem not set!"
        usage
        exit 1;
    fi
    if [ -z "$name" ];
    then
        echo "ERROR: -n/--name not set!"
        usage
        exit 1;
    fi
    if [ -z "$n_eigenvectors" ];
    then
        echo "ERROR: -z/--n_eigenvectors not set!"
        usage
        exit 1;
    fi
}

# [[ ${BASH_SOURCE[0]} = "$0" ]] -> main does not run when this script is sourced
# main "$@" -> Send the arguments to the main function (this way project flow can be at top)
# exit -> safeguard against the file being modified while it is interpreted
main "$@";






