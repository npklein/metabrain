#!/usr/bin/env bash

# This script clones the Molgenis Compute pipeline for alignin FastQ files and
# modifies the protocols and parameter files for analysis of Brain eQTL data.
# Makes samplesheet dependent parameter files.

# NOTE: This script changes paths in a parameters.csv that are specific to our clusters
#       If you want to use this somewhere else you have to change that in the script"


set -e
set -u

module load Python

cohort=
project_dir=
tmpdir=
results_dir=
script_dir=$(dirname "$0")

main(){
    parse_commandline "$@"
    # move to project dir if it exists.
    # I prefer making the project dir outside of this script in case wrong directory is given
    if [ -d "$project_dir" ];
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
    cohort_specific_steps
    make_pipeline_scripts
    cd -
}

usage(){
    # print the usage of the programme
    programname=$0
    echo "usage: $programname -c cohort -p project_directory -t tmp04"
    echo "  -c      provide the cohort to prepare pipeline for (mandatory)"
    echo "  -p      Base of the project_dir where jobs will be put and pipeline will be built (mandatory)"
    echo "  -r      Dir where results will be put (mandatory)"
    echo "  -t      tmpdir of the cluster (e.g. tmp03 for boxy, tmp04 for calculon). Will use /groups/$tmpdir/.../ (mandatory)"
    echo "  -h      display help"
    exit 1
}

clone_pipelines(){
    # Clone the relevant github repositories and move some files around
    echo "Cloning pipeline...."

    # In case the program has already run in current direcotry
    rm -rf Public_RNA-seq_QC
    rm -rf Public_RNA-seq_quantification
    rm -rf Public_RNA-seq_genotypeCalling
    git clone https://github.com/npklein/molgenis-pipelines.git
    # Move the relevant pipeline directories to the main project_dir
    mv molgenis-pipelines/compute5/Public_RNA-seq_QC/ .
    mv molgenis-pipelines/compute5/Public_RNA-seq_quantification/ .
    mv molgenis-pipelines/compute5/Public_RNA-seq_genotypeCalling/protocols/* Public_RNA-seq_QC/protocols/
    sort -u Public_RNA-seq_QC/parameter_files/parameters.csv molgenis-pipelines/compute5/Public_RNA-seq_genotypeCalling/parameter_files/parameters.csv > Public_RNA-seq_QC/parameter_files/parameters.csv.tmp
    mv Public_RNA-seq_QC/parameter_files/parameters.csv.tmp Public_RNA-seq_QC/parameter_files/parameters.csv
    mv molgenis-pipelines/compute5/Public_RNA-seq_genotypeCalling/chromosomes.csv Public_RNA-seq_QC/
    rsync -vP brain_eQTL/pipeline_prepare_scripts/modified_workflows/workflow_brain_eQTL.csv /tmp/workflow_brain_eQTL.csv
    rm -rf molgenis-pipelines
}

adjust_workflows(){
    # Adjust the molgenis compute workflows to only run the steps necesarry for this project and to add new steps
    echo "Adjusting workflows..."

    # Remove unused workflows to make workflow directory less messy
    rm Public_RNA-seq_QC/workflows/*
    mv /tmp/workflow_brain_eQTL.csv Public_RNA-seq_QC/workflows/

    # Remove all the steps we don't need from the workflow
    # At the moment there are none to remove

    # change some dependencies
    sed -i 's@HaplotypeCallerGvcf,../protocols/GatkHaplotypeCallerGvcf.sh,BQSR@HaplotypeCallerGvcf,../protocols/GatkHaplotypeCallerGvcf.sh,BQSR;AnalyseCovariates@' Public_RNA-seq_QC/workflows/workflow_brain_eQTL.csv
    sed -i 's@GatkSplitAndTrim,../protocols/GatkSplitAndTrim.sh,MarkDuplicates@GatkSplitAndTrim,../protocols/GatkSplitAndTrim.sh,MarkDuplicates;Flagstat;CollectMultipleMetrics;CollectRnaSeqMetrics@' Public_RNA-seq_QC/workflows/workflow_brain_eQTL.csv
    sed -i 's@MergeBamFiles,../protocols/MergeBamFiles.sh,AddReadGroup@MergeBamFiles,../protocols/MergeBamFiles.sh,AddReadGroup;CollectRnaSeqQcMetrics;CollectMultipleQcMetrics;GatkUnifiedGenotyper@' Public_RNA-seq_QC/workflows/workflow_brain_eQTL.csv

    # Change/add steps that are not in by default
    sed -i 's;Kallisto,../protocols/Kallisto.sh,;Sailfish,../protocols/Sailfish.sh,;' Public_RNA-seq_quantification/workflows/workflow.csv

    # add last step to remove bams
    echo "RemoveBamFiles,../protocols/RemoveBamFiles.sh,HaplotypeCallerGvcf" >> Public_RNA-seq_QC/workflows/workflow_brain_eQTL.csv
}



change_protocols(){
    # Adjust the molgenis compute protocols based on needs of the brain eQTL project
    echo "Changing protocols..."

    # Since we don't filter BAMs, change the input file in sortbam
    # has to be done before adding the removal of unfilteredbam to end of file so it doesnt replace
    sed -i 's;filteredBam;unfilteredBam;' Public_RNA-seq_QC/protocols/SortBam.sh

    # Although for ENA we do add in readgroup information as the new version of GATK needs it,
    # this was not done for the other cohorts. Therefore, remove  METRIC_ACCUMULATION_LEVEL here as well
    # Otherwise, tries to do it per readgroud. Now does it per file
    sed -i '/METRIC_ACCUMULATION_LEVEL/d' Public_RNA-seq_QC/protocols/CollectRnaSeqMetrics.sh

    # Original SAM to BAM conversion is done in seperate step, but this step is removed from the workflow
    # Add the conversion to the STAR alignment script
    echo "echo \"convert SAM to BAM\"" >> Public_RNA-seq_QC/protocols/STARMapping.sh
    echo "module load SAMtools/1.5-foss-2015b" >> Public_RNA-seq_QC/protocols/STARMapping.sh
    echo "mkdir -p ${results_dir}/results/\${unfilteredBamDir}/" >> Public_RNA-seq_QC/protocols/STARMapping.sh
    echo "samtools view -h -b \${alignmentDir}/\${uniqueID}.sam > \${unfilteredBamDir}/\${uniqueID}.bam" >> Public_RNA-seq_QC/protocols/STARMapping.sh
    echo "rm \${alignmentDir}/\${uniqueID}.sam" >> Public_RNA-seq_QC/protocols/STARMapping.sh

    # For quick genotyping to test relatedness we use GATK3, while later we want to use GATK4. So change in UnifiedGenotyper
    sed -i 's;GATK/${gatkVersion};GATK/3.8-0-Java-1.8.0_121;' Public_RNA-seq_QC/protocols/GatkUnifiedGenotyper.sh
    sed -i 's;GATK/${gatkVersion};GATK/3.8-0-Java-1.8.0_121;' Public_RNA-seq_QC/protocols/VariantEval.sh

    # MergeBam takes sortedBam instead of AddOrReplace Bams
    sed -i 's;addOrReplaceGroupsBam;sortedBam;' Public_RNA-seq_QC/protocols/MergeBamFiles.sh

    # add the removeBamFile protocol
    rsync -vP $script_dir/modified_protocols/RemoveBamFiles.sh Public_RNA-seq_QC/protocols/

    # All the GATK tools have been modified from the original to use GATK4, copy them over
    rsync -P $script_dir/modified_protocols/GatkAnalyseCovariates.sh Public_RNA-seq_QC/protocols/
    rsync -P $script_dir/modified_protocols/GatkBQSR.sh Public_RNA-seq_QC/protocols/
    rsync -P $script_dir/modified_protocols/GatkHaplotypeCallerGvcf.sh Public_RNA-seq_QC/protocols/
    rsync -P $script_dir/modified_protocols/GatkSplitAndTrim.sh Public_RNA-seq_QC/protocols/
}

change_parameter_files(){
    # Change those parameters in parameters.csv that are old (e.g. old path to resource files) or cluster specific (e.g. tmp03/tmp04)
    echo "Changing parameter files..."

    # Change the alignment pipeline parameter file
    # NOTE: The replacemets are cluster specific, but don't want to make this command line options because that would make too many of them
    #       Either change below code or change the parameters.csv file
    sed -i 's;group,umcg-wijmenga;group,umcg-biogen;' Public_RNA-seq_QC/parameter_files/parameters.csv
    sed -i "s;resDir,/groups/umcg-wijmenga/tmp04/resources/;resDir,/apps/data/;" Public_RNA-seq_QC/parameter_files/parameters.csv
    sed -i "s;resDir,/groups/umcg-wijmenga/tmp03/resources/;resDir,/apps/data/;" Public_RNA-seq_QC/parameter_files/parameters.csv
    sed -i "s;projectDir,\${root}/\${group}/\${tmp}/projects/umcg-ndeklein/\${project}/;projectDir,${results_dir}/results/;" Public_RNA-seq_QC/parameter_files/parameters.csv
    sed -i 's;STARindex;STARindex,${resDir}/ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/STAR/${starVersion}/;' Public_RNA-seq_QC/parameter_files/parameters.csv
    sed -i 's;alignmentDir,${projectDir}/hisat/;alignmentDir,${projectDir}/star;' Public_RNA-seq_QC/parameter_files/parameters.csv
    sed -i 's;fastqExtension,.gz;fastqExtension,.fq.gz;' Public_RNA-seq_QC/parameter_files/parameters.csv
    sed -i 's;goolf-1.7.20;foss-2015b;g' Public_RNA-seq_QC/parameter_files/parameters.csv
    sed -i 's;fastqExtension,.fq.gz;fastqExtension,.fastq.gz;' Public_RNA-seq_QC/parameter_files/parameters.csv
    sed -i 's;1.102-Java-1.7.0_80;1.119-Java-1.7.0_80;' Public_RNA-seq_QC/parameter_files/parameters.csv
    sed -i 's;onekgGenomeFasta,${resDir}/${genomeBuild}/indices/human_g1k_v${human_g1k_vers}.fasta;onekgGenomeFasta,${resDir}/ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/GRCh38.p5.genome.fa;' Public_RNA-seq_QC/parameter_files/parameters.csv
    sed -i 's;genesRefFlat,${resDir}/Ensembl/release-${ensemblVersion}/gtf/homo_sapiens/${genomeLatSpecies}.${genomeGrchBuild}.${ensemblVersion}.refflat;genesRefFlat,${resDir}/ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/gencode.v24.chr_patch_hapl_scaff.annotation.refflat;' Public_RNA-seq_QC/parameter_files/parameters.csv
    sed -i 's;rRnaIntervalList,${resDir}//picard-tools/Ensembl${ensemblVersion}/${genomeLatSpecies}.${genomeGrchBuild}.${ensemblVersion}.rrna.interval_list;rRnaIntervalList,${resDir}/ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/gencode.v24.chr_patch_hapl_scaff.annotation.rRNA.interval_list;' Public_RNA-seq_QC/parameter_files/parameters.csv
    sed -i 's;gatkVersion,3.4-0-Java-1.7.0_80;gatkVersion,4.0.8.1-foss-2015b-Python-3.6.3;' Public_RNA-seq_QC/parameter_files/parameters.csv
    sed -i 's;dbSNP/dbsnp_138.b37.vcf;ftp.ncbi.nlm.nih.gov/snp/organisms/archive/human_9606_b144_GRCh38p2/VCF/All_20150603.with_chr.vcf.gz;' Public_RNA-seq_QC/parameter_files/parameters.csv
    # Change in genotyping parameter file, but tht has been merged with the QC parameter file
    sed -i 's;splitAndTrimBam,${splitAndTrimDir}${sampleName}.bam;splitAndTrimBam,${splitAndTrimDir}${sampleName}.${chromosome}.bam;' Public_RNA-seq_QC/parameter_files/parameters.csv
    sed -i 's;bqsrAfterGrp,${bqsrDir}${sampleName}.after.grp;bqsrAfterGrp,${bqsrDir}${sampleName}.${chromosome}.after.grp;' Public_RNA-seq_QC/parameter_files/parameters.csv
    sed -i 's;bqsrBeforeGrp,${bqsrDir}${sampleName}.before.grp;bqsrBeforeGrp,${bqsrDir}${sampleName}.${chromosome}.before.grp;' Public_RNA-seq_QC/parameter_files/parameters.csv
    sed -i 's;bqsrBam,${bqsrDir}${sampleName}.bam;bqsrBam,${bqsrDir}${sampleName}.${chromosome}.bam;' Public_RNA-seq_QC/parameter_files/parameters.csv
    sed -i 's;${genomeBuild}/sv/1000G/Mills_and_1000G_gold_standard.indels.b37.vcf;storage.cloud.google.com/genomics-public-data/resources/broard/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz;' Public_RNA-seq_QC/parameter_files/parameters.csv
    sed -i 's;${genomeBuild}/sv/1000G/1000G_phase1.indels.b37.vcf;storage.cloud.google.com/genomics-public-data/resources/broard/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz;' Public_RNA-seq_QC/parameter_files/parameters.csv

    # Change the qunatification pipeline parameter file
    sed -i 's;group,umcg-wijmenga;group,umcg-biogen;' Public_RNA-seq_quantification/parameter_files/parameters.csv
    sed -i "s;tmp03;$tmpdir;" Public_RNA-seq_quantification/parameter_files/parameters.csv
    sed -i "s;projectDir,\${root}/\${group}/\${tmp}/projects/umcg-ndeklein/\${project}/;projectDir,${results_dir}/results/;" Public_RNA-seq_quantification/parameter_files/parameters.csv
    sed -i "s;projects/umcg-ndeklein/\${project};biogen/input/${cohort}/results/;" Public_RNA-seq_quantification/parameter_files/parameters.csv
    chr38gtf="/apps/data/ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/gencode.v24.chr_patch_hapl_scaff.annotation.gtf.gz"
    sed -i "s;annotationGtf,/apps/data/ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf;annotationGtf,${chr38gtf};" Public_RNA-seq_quantification/parameter_files/parameters.csv
    echo "bedtoolsVersion,2.25.0-foss-2015b" >> Public_RNA-seq_quantification/parameter_files/parameters.csv
    echo "picardVersion,2.10.0-foss-2015b-Java-1.8.0_74" >> Public_RNA-seq_quantification/parameter_files/parameters.csv
    echo "iolibVersion,1.14.6-foss-2015b" >> Public_RNA-seq_quantification/parameter_files/parameters.csv
    echo "resDir,/apps/data/" >> Public_RNA-seq_quantification/parameter_files/parameters.csv
    echo 'onekgGenomeFasta,${resDir}/ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/GRCh38.p5.genome.fa' >> Public_RNA-seq_quantification/parameter_files/parameters.csv
    sed -i 's;htseqVersion,0.6.1p1-foss-2015b;htseqVersion,0.9.1-foss-2015b-Python-2.7.11;' Public_RNA-seq_quantification/parameter_files/parameters.csv
    # most of the protocols are not stranded, so change here. If protocol is stranded for one of the
    # cohorts it will be adjusted in cohort specific method below
    sed -i 's;stranded,reverse;stranded,no;' Public_RNA-seq_quantification/parameter_files/parameters.csv


    # Change in both
    sed -i 's;genomeBuild,b37;genomeBuild,b38;' Public_RNA-seq*/parameter_files/parameters.csv
    sed -i 's;genomeGrchBuild,GRCh37;genomeGrchBuild,GRCh38;' Public_RNA-seq*/parameter_files/parameters.csv
    sed -i 's;human_g1k_vers,37;human_g1k_vers,38;' Public_RNA-seq*/parameter_files/parameters.csv
    sed -i 's;ensemblVersion,75;ensemblVersion,?;' Public_RNA-seq*/parameter_files/parameters.csv
    sed -i 's;${resDir}/UMCG/1000G_interval_list//ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.MAF_0.05.SNPs_only.recode.annotated.EXON_only.interval_list;${resDir}/ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/gencode.v24.chr_patch_hapl_scaff.annotation.exon.interval_list;' Public_RNA-seq*/parameter_files/parameters.csv
}


make_samplesheets(){
    # Make the samplesheets. How the samplesheet is made is dependent om the cohort
    echo "Making samplesheets..."
    rm Public_RNA-seq_QC/samplesheet1.csv
    mkdir -p Public_RNA-seq_QC/samplesheets/

    samplesheet_script_dir=$script_dir/samplesheet_scripts/
    if [[ "$cohort" == "TargetALS" ]];
    then
        echo "ERROR: Don't need to get genotypes for the TargetALS samples. Because the pipeline needs to be set up quite differently, use make_alignmentAndQuantification_pipeline.sh instead"
        exit 1;
    elif [[ "$cohort" == "CMC" ]];
    then
        echo "ERROR: Don't need to get genotypes for the CMC samples. Because the pipeline needs to be set up quite differently, use make_alignmentAndQuantification_pipeline.sh instead"
        exit 1;
    elif [[ "$cohort" == "Braineac" ]];
    then
        echo "ERROR: Don't need to get genotypes for the Braineac samples. Because the pipeline needs to be set up quite differently, use make_alignmentAndQuantification_pipeline.sh instead"
        exit 1;
    elif [[ "$cohort" == "ENA" ]];
    then
        # Need to know where to find samplesheet_ENA_20181212.tx, this is in the brain_eQTL github directory. Since this script is also in this directory, find the directory like so (and add it to parameters file):
        brain_eQTL_dir="$( cd "$( dirname $( dirname "${BASH_SOURCE[0]}" ) )" >/dev/null && pwd )"
        python $samplesheet_script_dir/make_samplesheet_ENA.py $brain_eQTL_dir/ENA/ENA_samplesheets/samplesheet_ENA_20181212.txt \
                                                               /scratch/umcg-ndeklein/tmp03/ENA/pipelines/results/fastq/
    else
        echo "No code written for cohort $cohort"
        exit 1;
    fi
}

change_prepare_scripts(){
    echo "Changing prepare scripts..."

    # Do general changes for the alignment and haplotyping pipeline
    sed -i 's;/path/to/molgenis-pipelines/compute5/Public_RNA-seq_QC/;;' Public_RNA-seq_QC/prepare_Public_RNA-seq_QC.sh
    sed -i "s;/groups/umcg-wijmenga/tmp04/umcg-ndeklein/molgenis-pipelines/compute5/Public_RNA-seq_QC/;;" Public_RNA-seq_QC/prepare_Public_RNA-seq_QC.sh
    sed -i "s;parameters.converted.csv;${project_dir}/Public_RNA-seq_QC/parameter_files/parameters.converted.csv;" Public_RNA-seq_QC/prepare_Public_RNA-seq_QC.sh
    sed -i "s;workflow.csv;${project_dir}/Public_RNA-seq_QC/workflows//workflow_brain_eQTL.csv;" Public_RNA-seq_QC/prepare_Public_RNA-seq_QC.sh
    sed -i "s;-rundir /groups/umcg-wijmenga/tmp04/umcg-ndeklein/rundirs/QC/;-rundir ${project_dir}/jobs/alignmentAndHaplotyping/alignmentDir;" Public_RNA-seq_QC/prepare_Public_RNA-seq_QC.sh
    sed -i "s;/groups/umcg-wijmenga/tmp04/umcg-ndeklein/samplesheets/;${project_dir}/Public_RNA-seq_QC/samplesheets/;" Public_RNA-seq_QC/prepare_Public_RNA-seq_QC.sh


    sed -i "s;--generate \\\;--generate -p ${project_dir}/Public_RNA-seq_QC/chromosomes.csv \\\;" Public_RNA-seq_QC/prepare_Public_RNA-seq_QC.sh
    # add chr in front of the chromosomes.csv file to follow b38 rules
    awk '{print "chr" $0;}' Public_RNA-seq_QC/chromosomes.csv > Public_RNA-seq_QC/chromosomes.csv.tmp
    sed 's;chrchromosome;chromosome;' Public_RNA-seq_QC/chromosomes.csv.tmp > Public_RNA-seq_QC/chromosomes.csv
    rm Public_RNA-seq_QC/chromosomes.csv.tmp
    # replace the chr23, 24 and 25 for chrX, chrY, and chrM
    sed -i 's;chr23;chrX;' Public_RNA-seq_QC/chromosomes.csv
    sed -i 's;chr24;chrY;' Public_RNA-seq_QC/chromosomes.csv
    sed -i 's;chr25;chrM;' Public_RNA-seq_QC/chromosomes.csv

    # Do general changes for the quantification pipeline
    sed -i "s;workflows/workflow.csv;${project_dir}/Public_RNA-seq_quantification/workflows/workflow.csv;" Public_RNA-seq_quantification/prepare_quantification.sh
    sed -i "s;parameter_files/;${project_dir}/Public_RNA-seq_quantification/parameter_files/;" Public_RNA-seq_quantification/prepare_quantification.sh
    sed -i "s;-rundir ${project_dir}/jobs/alignment/alignmentDir;-rundir ${project_dir}/jobs/quantification/alignmentDir;" Public_RNA-seq_QC/prepare_Public_RNA-seq_QC.sh


    # Do specific changes per samplesheet batch
    mkdir Public_RNA-seq_QC/prepare_scripts/
    mkdir Public_RNA-seq_quantification/prepare_scripts/
    for samplesheet in `ls Public_RNA-seq_QC/samplesheets/samplesheet*txt`;
    do
        name=$(echo ${samplesheet} | awk -F"." '{print $2}')

        # Change it for the alignment pipeline per batch
        echo "making prepare script using $samplesheet: Public_RNA-seq_QC/prepare_scripts/prepare_Public_RNA-seq_QC.$name.sh"
        rsync -P Public_RNA-seq_QC/prepare_Public_RNA-seq_QC.sh Public_RNA-seq_QC/prepare_scripts/prepare_Public_RNA-seq_QC.$name.sh
        sed -i "s;samplesheet.csv;samplesheet_${cohort}_RNA.$name.txt;" Public_RNA-seq_QC/prepare_scripts/prepare_Public_RNA-seq_QC.$name.sh
        sed -i "s;alignmentDir;$name;" Public_RNA-seq_QC/prepare_scripts/prepare_Public_RNA-seq_QC.$name.sh

        # Change it for the quantificaiton pipeline per batch
        rsync -P Public_RNA-seq_quantification/prepare_quantification.sh Public_RNA-seq_quantification/prepare_scripts/prepare_quantification.$name.sh
        sed -i "s;samplesheet.csv;${project_dir}/Public_RNA-seq_QC/samplesheets/samplesheet_${cohort}_RNA.$name.txt;" Public_RNA-seq_quantification/prepare_scripts/prepare_quantification.$name.sh
        sed -i "s;alignmentDir;$name;" Public_RNA-seq_quantification/prepare_scripts/prepare_quantification.$name.sh
        break
    done
}


make_pipeline_scripts(){
    # Make the pipeline scripts using Molgenis Compute
    echo "Maing pipeline scripts..."
    echo $PWD

    # convert the parameter files to wide format now instead of in the change_parameters function because the cohort_specific_steps function might have changed some of them
    bash Public_RNA-seq_QC/parameter_files/convert.sh Public_RNA-seq_QC/parameter_files/parameters.csv Public_RNA-seq_QC/parameter_files/parameters.converted.csv
    bash Public_RNA-seq_QC/parameter_files/convert.sh Public_RNA-seq_quantification/parameter_files/parameters.csv Public_RNA-seq_quantification/parameter_files/parameters.converted.csv
    cd Public_RNA-seq_QC/prepare_scripts;
    for f in *sh;
    do
        bash $f;
    done
    cd -


#    cd Public_RNA-seq_quantification/prepare_scripts
#    for f in *sh;
#    do
#        bash $f;
#    done
#    cd -
}


cohort_specific_steps(){
    # Collection of cohort specific steps
    echo "Running cohort specific steps..."

    if [[ "$cohort" == "ENA" ]];
    then
        echo "ENA specific methods..."
        # Need to know where to find download_ENA_samples.py, this is in the brain_eQTL github directory. Since this script is also in this directory, find the directory like so (and add it to parameters file):
        brain_eQTL_dir="$( cd "$( dirname $( dirname "${BASH_SOURCE[0]}" ) )" >/dev/null && pwd )"
        echo "brainDir,$brain_eQTL_dir" >> Public_RNA-seq_QC/parameter_files/parameters.csv
        echo "enaSamplesheet,$brain_eQTL_dir/ENA/ENA_samplesheets/samplesheet_ENA_20181212.txt" >> Public_RNA-seq_QC/parameter_files/parameters.csv
        # Because we only need to download from ENA if the cohort is ENA, copy over
        rsync -P $script_dir/modified_protocols/DownloadFromENA.sh Public_RNA-seq_QC/protocols/

        # Different molgenis compute script is installed on Peregrine where the public data is run, so change
        for prepare_script in Public_RNA-seq*/prepare_scripts/*sh;
        do
            sed -i 's;v16.05.1-Java-1.8.0_45;v16.11.1-Java-1.8.0_74;' $prepare_script
        done

        # have to change a lot of parameters because the script locations are all different
        sed -i 's;WORKDIR,/groups/;WORKDIR,/scratch/;' Public_RNA-seq*/parameter_files/parameters.csv
        sed -i 's;group,umcg-biogen;group,umcg-ndeklein;' Public_RNA-seq*/parameter_files/parameters.csv
        sed -i "s;resDir,/apps/data/;resDir,/data/umcg-ndeklein/apps/data/;" Public_RNA-seq_QC/parameter_files/parameters.csv
    fi
}

parse_commandline(){
    # Parse the command line arguments
    echo "Parsing commandline..."

    # Check to see if at least one argument is given
    if [ $# -eq 0 ]
    then
        echo "ERROR: No arguments supplied"
        usage
        exit 1;
    fi

    while [[ $# -ge 1 ]]; do
        case $1 in
            -c | --cohort )         shift
                                    cohort=$1
                                    ;;
            -p | --project_dir )        shift
                                    project_dir=$1
                                    ;;
            -r | --results_dir )    shift
                                    results_dir=$1
                                    ;;
            -t | --tmpdir )         shift
                                    tmpdir=$1
                                    ;;
            -h | --help )           usage
                                    exit
                                    ;;
            * )                     echo "ERROR: Undexpected argument: $1"
                                    usage
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
    if [ -z "$tmpdir" ];
    then
        echo "ERROR: -t/--tmpdir not set!"
        usagae
        exit 1;
    fi
    if [ -z "$results_dir" ];
    then
        echo "ERROR: -r/--results_dir not set!"
        usagae
        exit 1;
    fi
    echo "Making pipeline..."
}

# [[ ${BASH_SOURCE[0]} = "$0" ]] -> main does not run when this script is sourced
# main "$@" -> Send the arguments to the main function (this way project flow can be at top)
# exit -> safeguard against the file being modified while it is interpreted
[[ ${BASH_SOURCE[0]} = "$0" ]] && main "$@"; exit;
