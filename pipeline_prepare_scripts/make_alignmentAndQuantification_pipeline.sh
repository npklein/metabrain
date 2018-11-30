#!/usr/bin/env bash

# Download and make quantification pipeline
# This script clones the Molgenis Compute pipeline for alignin FastQ files and
# modifies the protocols and parameter files for analysis of Brain eQTL data.
# Makes samplesheet dependent parameter files.
# Author: Niek de Klein

set -e
set -u

module load Python
cohort=
project_dir=

main(){
    parse_commandline "$@"
    # move to project dir if it exists.
    # I prefer making the project dir outside of this script in case wrong directory is given
    if [ -d "$DIRECTORY" ];
    then
        echo "Changing directory to $project_dir"
        cd $project_dir
    else
        echo "ERROR: $project_dir does not exist"
        exit 1;
    fi
    clone_pipelines
    adjust_workflows
    change_protocols
    change_parameter_files
    make_samplesheets
    change_prepare_scripts
    make_pipeline_scripts
    cd -
}

usage(){
    # print the usage of the programme
    programname=$0
    echo "Download and make qunatification pipeline."
    echo "NOTE: This script changes paths in a parameters.csv that are specific to our clusters"
    echo "      If you want to use this somewhere else you have to change that in the script"
    echo "usage: $programname -c cohort -p project_directory"
    echo "  -c      provide the cohort to prepare pipeline for (mandatory)"
    echo "  -p      Base of the project_dir where pipeline scripts will be put (mandatory)"
    echo "  -h      display help"
    exit 1
}

clone_pipelines(){
    # Clone the relevant github repositories and move some files around
    echo "Cloning pipeline...."

    # In case the program has already run in current direcotry
    rm -rf Public_RNA-seq_QC
    rm -rf Public_RNA-seq_quantification
    git clone https://github.com/npklein/molgenis-pipelines.git

    # Move the relevant pipeline directories to the main project_dir
    mv molgenis-pipelines/compute5/Public_RNA-seq_QC/ .
    mv molgenis-pipelines/compute5/Public_RNA-seq_quantification/ .
    rm -rf molgenis-pipelines
}

adjust_workflows(){
    # Adjust the molgenis compute workflows to only run the steps necesarry for this project and to add new steps
    echo "Adjusting workflows..."

    # Remove unused workflows to make workflow directory less messy
    mv Public_RNA-seq_QC/workflows/workflowSTAR.csv /tmp/workflowSTAR.csv
    rm Public_RNA-seq_QC/workflows/*
    mv /tmp/workflowSTAR.csv Public_RNA-seq_QC/workflows/

    # Remove all the steps we don't need from the workflow
    sed -i '/VerifyBamID/d' Public_RNA-seq_QC/workflows/workflowSTAR.csv
    sed -i '/VariantEval/d' Public_RNA-seq_QC/workflows/workflowSTAR.csv
    sed -i '/HtseqCount/d' Public_RNA-seq_quantification/workflows/workflow.csv
    sed -i '/GatkUnifiedGenotyper/d' Public_RNA-seq_QC/workflows/workflowSTAR.csv
    sed -i '/SortBam.sh/d' Public_RNA-seq_QC/workflows/workflowSTAR.csv
    sed -i '/AddOr/d' Public_RNA-seq_QC/workflows/workflowSTAR.csv

    # Change/add steps that are not in by default
    sed -i 's/SortBam/CreateCramFiles/' Public_RNA-seq_QC/workflows/workflowSTAR.csv
    sed -i 's;Kallisto,../protocols/Kallisto.sh,;Sailfish,../protocols/Sailfish.sh,;' Public_RNA-seq_quantification/workflows/workflow.csv
    echo "FeatureCount,../protocols/FeatureCount.sh," >> Public_RNA-seq_quantification/workflows/workflow.csv
}



change_protocols(){
    # Adjust the molgenis compute protocols based on needs of the brain eQTL project
    echo "Changing protocols..."

    # add a line to delete the unfiltered BAM after cramming
    echo "echo \"remove \${unfilteredBamDir}/\${uniqueID}.bam\"" >> Public_RNA-seq_QC/protocols/CreateCramFiles.sh
    echo "rm \${unfilteredBamDir}/\${uniqueID}.bam" >> Public_RNA-seq_QC/protocols/CreateCramFiles.sh

    # Because we first convert to cram before running the collectMetrics jobs, change this for all Collect*sh scripts
    sed -i 's;#string sortedBam;#string cramFileDir;' Public_RNA-seq_QC/protocols/Collect*sh
    sed -i 's;#string sortedBai;#string uniqueID;' Public_RNA-seq_QC/protocols/Collect*sh
    sed -i 's;${sortedBam};${cramFileDir}${uniqueID}.cram;' Public_RNA-seq_QC/protocols/Collect*sh

    # Did not add readgroup information at this point, so remove line METRIC_ACCUMULATION_LEVEL
    # Otherwise, tries to do it per readgroud. Now does it per file
    sed -i '/METRIC_ACCUMULATION_LEVEL/d' Public_RNA-seq_QC/protocols/CollectRnaSeqMetrics.sh

    # Original SAM to BAM conversion is done in seperate step, but this step is removed from the workflow
    # Add the conversion to the STAR alignment script
    echo "echo \"convert SAM to BAM\"" >> Public_RNA-seq_QC/protocols/STARMapping.sh
    echo "module load SAMtools/1.5-foss-2015b" >> Public_RNA-seq_QC/protocols/STARMapping.sh
    echo "mkdir -p /groups/umcg-biogen/tmp04/biogen/input/TargetALS/pipelines/results/\${unfilteredBamDir}/" >> Public_RNA-seq_QC/protocols/STARMapping.sh
    echo "samtools view -h -b \${alignmentDir}/\${uniqueID}.sam > \${unfilteredBamDir}/\${uniqueID}.bam" >> Public_RNA-seq_QC/protocols/STARMapping.sh
    echo "rm \${alignmentDir}/\${uniqueID}.sam" >> Public_RNA-seq_QC/protocols/STARMapping.sh

    # Removed the filteredBam step, so change this in the SortBam protocol
    sed -i 's;filteredBam;unfilteredBam;' Public_RNA-seq_QC/protocols/SortBam.sh
}

change_parameter_files(){
    # Change those parameters in parameters.csv that are old (e.g. old path to resource files) or cluster specific (e.g. tmp03/tmp04)
    echo "Changing parameter files..."

    # Change the alignment pipeline parameter file
    # NOTE: The replacemets are cluster specific, but don't want to make this command line options because that would make too many of them
    #       Either change below code or change the parameters.csv file
    sed -i 's;group,umcg-wijmenga;group,umcg-biogen;' Public_RNA-seq_QC/parameter_files/parameters.csv
    sed -i 's;resDir,/groups/umcg-wijmenga/tmp04/resources/;resDir,/apps/data/;' Public_RNA-seq_QC/parameter_files/parameters.csv
    sed -i 's;projects/umcg-ndeklein/${project};biogen/input/TargetALS/pipelines/results/;' Public_RNA-seq_QC/parameter_files/parameters.csv
    sed -i 's;STARindex;STARindex,${resDir}/ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/STAR/${starVersion}/;' Public_RNA-seq_QC/parameter_files/parameters.csv
    sed -i 's;alignmentDir,${project_dir}/hisat/;alignmentDir,${project_dir}/star;' Public_RNA-seq_QC/parameter_files/parameters.csv
    sed -i 's;fastqExtension,.gz;fastqExtension,.fq.gz;' Public_RNA-seq_QC/parameter_files/parameters.csv
    sed -i 's;goolf-1.7.20;foss-2015b;g' Public_RNA-seq_QC/parameter_files/parameters.csv
    sed -i 's;fastqExtension,.fq.gz;fastqExtension,.fastq.gz;' Public_RNA-seq_QC/parameter_files/parameters.csv
    sed -i 's;1.102-Java-1.7.0_80;1.119-Java-1.7.0_80;' Public_RNA-seq_QC/parameter_files/parameters.csv
    sed -i 's;onekgGenomeFasta,${resDir}/${genomeBuild}/indices/human_g1k_v${human_g1k_vers}.fasta;onekgGenomeFasta,${resDir}/ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/GRCh38.p5.genome.fa;' Public_RNA-seq_QC/parameter_files/parameters.csv
    sed -i 's;genesRefFlat,${resDir}/picard-tools/Ensembl${ensemblVersion}/${genomeLatSpecies}.${genomeGrchBuild}.${ensemblVersion}.refflat;genesRefFlat,${resDir}/ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/gencode.v24.chr_patch_hapl_scaff.annotation.refflat;' Public_RNA-seq_QC/parameter_files/parameters.csv
    sed -i 's;rRnaIntervalList,${resDir}//picard-tools/Ensembl${ensemblVersion}/${genomeLatSpecies}.${genomeGrchBuild}.${ensemblVersion}.rrna.interval_list;rRnaIntervalList,${resDir}/ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/gencode.v24.chr_patch_hapl_scaff.annotation.rRNA.interval_list;' Public_RNA-seq_QC/parameter_files/parameters.csv
    bash Public_RNA-seq_QC/parameter_files/convert.sh Public_RNA-seq_QC/parameter_files/parameters.csv Public_RNA-seq_QC/parameter_files/parameters.converted.csv

    # Change the qunatification pipeline parameter file
    sed -i 's;group,umcg-wijmenga;group,umcg-biogen;' Public_RNA-seq_quantification/parameter_files/parameters.csv
    sed -i 's;tmp03;tmp04;' Public_RNA-seq_quantification/parameter_files/parameters.csv
    sed -i 's;projects/umcg-ndeklein/${project};biogen/input/TargetALS/pipelines/results/;' Public_RNA-seq_quantification/parameter_files/parameters.csv
    bash Public_RNA-seq_QC/parameter_files/convert.sh Public_RNA-seq_quantification/parameter_files/parameters.csv Public_RNA-seq_quantification/parameter_files/parameters.converted.csv
}


make_samplesheets(){
    # Make the samplesheets. How the samplesheet is made is dependent om the cohort
    echo "Making samplesheets..."

    rm Public_RNA-seq_QC/samplesheet1.csv
    mkdir -p Public_RNA-seq_QC/samplesheets/
    python make_samplesheet.py sample_annotation/UMG-Target_ALS_RNA_Clinical_Data_06122018.txt Public_RNA-seq_QC/samplesheets/samplesheet_TargetALS_RNA.
    python make_samplesheet.py sample_annotation/RNA_Metadata_TALS_2_5July2018.txt Public_RNA-seq_QC/samplesheets/samplesheet_TargetALS_RNA.samplesheet5july2018_
}

change_prepare_scripts(){
    # Change the molgenis compute prepare scirpts
    echo "Chaning prepare scripts..."

    # Do general changes for the alignment pipeline
    sed -i 's;/path/to/molgenis-pipelines/compute5/Public_RNA-seq_QC/;;' Public_RNA-seq_QC/prepare_Public_RNA-seq_QC.sh
    sed -i 's;/groups/umcg-wijmenga/tmp04/umcg-ndeklein/molgenis-pipelines/compute5/Public_RNA-seq_QC/;;' Public_RNA-seq_QC/prepare_Public_RNA-seq_QC.sh
    sed -i 's;parameters.converted.csv;/groups/umcg-biogen/tmp04/biogen/input/TargetALS/Public_RNA-seq_QC/parameter_files/parameters.converted.csv;' Public_RNA-seq_QC/prepare_Public_RNA-seq_QC.sh
    sed -i 's;workflow.csv;/groups/umcg-biogen/tmp04/biogen/input/TargetALS/Public_RNA-seq_QC/workflows/workflowSTAR.csv;' Public_RNA-seq_QC/prepare_Public_RNA-seq_QC.sh
    sed -i 's;-rundir /groups/umcg-wijmenga/tmp04/umcg-ndeklein/rundirs/QC/;-rundir /groups/umcg-biogen/tmp04/biogen/input/TargetALS/pipelines/alignment/alignmentDir;' Public_RNA-seq_QC/prepare_Public_RNA-seq_QC.sh
    sed -i 's;/groups/umcg-wijmenga/tmp04/umcg-ndeklein/samplesheets/;/groups/umcg-biogen/tmp04/biogen/input/TargetALS/Public_RNA-seq_QC/samplesheets/;' Public_RNA-seq_QC/prepare_Public_RNA-seq_QC.sh

    # Do general changes for the quantification pipeline
    sed -i 's;workflows/workflow.csv;/groups/umcg-biogen/tmp04/biogen/input/TargetALS/Public_RNA-seq_quantification/workflows/workflow.csv;' Public_RNA-seq_quantification/prepare_quantification.sh
    sed -i 's;parameter_files/;/groups/umcg-biogen/tmp04/biogen/input/TargetALS/Public_RNA-seq_quantification/parameter_files/;' Public_RNA-seq_quantification/prepare_quantification.sh


    # Do specific changes per samplesheet batch
    mkdir Public_RNA-seq_QC/prepare_scripts/
    mkdir Public_RNA-seq_quantification/prepare_scripts/
    for samplesheet in `ls Public_RNA-seq_QC/samplesheets/samplesheet*txt`;
    do
        name=$(echo ${samplesheet} | awk -F"." '{print $2}')

        # Change it for the alignment pipeline per batch
        echo "making prepare script using $samplesheet: Public_RNA-seq_QC/prepare_scripts/prepare_Public_RNA-seq_QC.$name.sh"
        rsync -vP Public_RNA-seq_QC/prepare_Public_RNA-seq_QC.sh Public_RNA-seq_QC/prepare_scripts/prepare_Public_RNA-seq_QC.$name.sh
        sed -i "s;samplesheet.csv;samplesheet_TargetALS_RNA.$name.txt;" Public_RNA-seq_QC/prepare_scripts/prepare_Public_RNA-seq_QC.$name.sh
        sed -i "s;alignmentDir;$name;" Public_RNA-seq_QC/prepare_scripts/prepare_Public_RNA-seq_QC.$name.sh

        # Change it for the quantificaiton pipeline per batch
        rsync -vP Public_RNA-seq_quantification/prepare_quantification.sh Public_RNA-seq_quantification/prepare_scripts/prepare_quantification.$name.sh
        sed -i "s;samplesheet.csv;/groups/umcg-biogen/tmp04/biogen/input/TargetALS/Public_RNA-seq_QC/samplesheets/samplesheet_TargetALS_RNA.$name.txt;" Public_RNA-seq_quantification/prepare_scripts/prepare_quantification.$name.sh
        sed -i "s;alignmentDir;$name;" Public_RNA-seq_quantification/prepare_scripts/prepare_quantification.$name.sh
        sed -i "s;-rundir results/;-rundir /groups/umcg-biogen/tmp04/biogen/input/TargetALS/pipelines/quantification/$name;" Public_RNA-seq_quantification/prepare_scripts/prepare_quantification.$name.sh
    done
}


make_pipeline_scripts(){
    # Make the pipeline scripts using Molgenis Compute
    echo "Maing pipeline scripts..."

    cd Public_RNA-seq_QC/prepare_scripts;
    for f in *sh;
    do
        bash $f;
    done
    cd -

    cd Public_RNA-seq_quantification/prepare_scripts
    for f in *sh;
    do
        bash $f;
    done
    cd -
}


cohort_specific_steps(){
    # Collection of cohort specific steps
    echo "Running cohort specific steps..."

#    rsync -vP STARMappingTwoPass.sh Public_RNA-seq_QC/protocols/
#   echo "STARMappingTwoPass,../protocols/STARMappingTwoPass.sh,CreateCramFiles" >> Public_RNA-seq_QC/workflows/workflowSTAR.csv

}

parse_commandLine(){
    echo "Parsing commandline..."

    while [ "$1" != "" ]; do
        case $1 in
            -c | --cohort )         shift
                                    cohort=$1
                                    ;;
            -p | --project_dir )     shift
                                    project_dir=$1
                                    ;;
            -h | --help )           usage
                                    exit
                                    ;;
            * )                     usage
                                    exit 1
        esac
        shift
    done
    echo "Cohort = $cohort"
    echo "Project directory = $project_dir"

    # if -z tests if variable is empty. Make sure the relevant variables are set
    if [ -z "$cohort" ];
    then
        echo "ERROR: -c/--cohort not set!"
        usage
        exit 1;
    fi
    if [ -z "$project_dir" ];
    then
        echo "ERROR: -p/--project_dir not set!"
        usagae
        exit 1;
    fi

    echo "Making pipeline..."
}


main "$@"
