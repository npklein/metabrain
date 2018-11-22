set -e
set -u

if [[ $# -eq 0 ]] ;
then
    echo 'No arguments given, cohort name has to be first argument'
    exit 1
fi

module load Python

rm -rf Public_RNA-seq_QC
rm -rf Public_RNA-seq_quantification

echo "Making pipeline..."
echo "Clone from github..."
git clone https://github.com/npklein/molgenis-pipelines.git

echo "Adjust workflow for out preference..."
mv molgenis-pipelines/compute5/Public_RNA-seq_QC/ .
mv molgenis-pipelines/compute5/Public_RNA-seq_quantification/ .
rm -rf molgenis-pipelines
mv Public_RNA-seq_QC/workflows/workflowSTAR.csv /tmp/workflowSTAR.csv
rm Public_RNA-seq_QC/workflows/*
mv /tmp/workflowSTAR.csv Public_RNA-seq_QC/workflows/
rsync -vP STARMappingTwoPass.sh Public_RNA-seq_QC/protocols/
sed -i '/VerifyBamID/d' Public_RNA-seq_QC/workflows/workflowSTAR.csv
sed -i '/VariantEval/d' Public_RNA-seq_QC/workflows/workflowSTAR.csv
sed -i '/HtseqCount/d' Public_RNA-seq_quantification/workflows/workflow.csv
sed -i '/GatkUnifiedGenotyper/d' Public_RNA-seq_QC/workflows/workflowSTAR.csv
sed -i '/SortBam.sh/d' Public_RNA-seq_QC/workflows/workflowSTAR.csv
sed -i '/AddOr/d' Public_RNA-seq_QC/workflows/workflowSTAR.csv
sed -i 's/SortBam/CreateCramFiles/' Public_RNA-seq_QC/workflows/workflowSTAR.csv
echo "STARMappingTwoPass,../protocols/STARMappingTwoPass.sh,CreateCramFiles" >> Public_RNA-seq_QC/workflows/workflowSTAR.csv
sed -i 's;Kallisto,../protocols/Kallisto.sh,;Sailfish,../protocols/Sailfish.sh,;' Public_RNA-seq_quantification/workflows/workflow.csv
echo "DONE"

echo "Changing protocols"
echo "delete dbsnp"
sed -i '/--dbsnp/d' Public_RNA-seq_QC/protocols/GatkUnifiedGenotyper.sh
echo "add a line to delete the unfiltered BAM after cramming"
echo "echo \"remove \${unfilteredBamDir}/\${uniqueID}.bam\"" >> Public_RNA-seq_QC/protocols/CreateCramFiles.sh
echo "rm \${unfilteredBamDir}/\${uniqueID}.bam" >> Public_RNA-seq_QC/protocols/CreateCramFiles.sh
sed -i 's;#string sortedBam;#string cramFileDir;' Public_RNA-seq_QC/protocols/Collect*sh
sed -i 's;#string sortedBai;#string uniqueID;' Public_RNA-seq_QC/protocols/Collect*sh
sed -i 's;${sortedBam};${cramFileDir}${uniqueID}.cram;' Public_RNA-seq_QC/protocols/Collect*sh
sed -i '/METRIC_ACCUMULATION_LEVEL/d' Public_RNA-seq_QC/protocols/CollectRnaSeqMetrics.sh
echo "add a line to delete the sorted BAM after unified genotyper"
echo """echo \"remove all sorted bam files\"
for sortedBamFile in \"\${sortedBam[@]}\"
do
    echo \"remove \$sortedBamFile\"
    rm \"\$sortedBamFile\"
done""" >> Public_RNA-seq_QC/protocols/VerifyBamID.sh
echo "DONE"
echo "add sam to bam to star"
echo "echo \"convert SAM to BAM\"" >> Public_RNA-seq_QC/protocols/STARMapping.sh
echo "module load SAMtools/1.5-foss-2015b" >> Public_RNA-seq_QC/protocols/STARMapping.sh
echo "mkdir -p /groups/umcg-biogen/tmp04/biogen/input/$1/pipelines/results/\${unfilteredBamDir}/" >> Public_RNA-seq_QC/protocols/STARMapping.sh
echo "samtools view -h -b \${alignmentDir}/\${uniqueID}.sam > \${unfilteredBamDir}/\${uniqueID}.bam" >> Public_RNA-seq_QC/protocols/STARMapping.sh
echo "rm \${alignmentDir}/\${uniqueID}.sam" >> Public_RNA-seq_QC/protocols/STARMapping.sh
echo "change filterBam to unfilteredBam in SortBam"
sed -i 's;filteredBam;unfilteredBam;' Public_RNA-seq_QC/protocols/SortBam.sh


echo "Changing parameters file"
sed -i 's;group,umcg-wijmenga;group,umcg-biogen;' Public_RNA-seq_QC/parameter_files/parameters.csv
sed -i 's;resDir,/groups/umcg-wijmenga/tmp04/resources/;resDir,/apps/data/;' Public_RNA-seq_QC/parameter_files/parameters.csv
sed -i 's;projects/umcg-ndeklein/${project};biogen/input/$1/pipelines/results/;' Public_RNA-seq_QC/parameter_files/parameters.csv
sed -i 's;STARindex;STARindex,${resDir}/ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/STAR/${starVersion}/;' Public_RNA-seq_QC/parameter_files/parameters.csv
sed -i 's;alignmentDir,${projectDir}/hisat/;alignmentDir,${projectDir}/star;' Public_RNA-seq_QC/parameter_files/parameters.csv
sed -i 's;fastqExtension,.gz;fastqExtension,.fq.gz;' Public_RNA-seq_QC/parameter_files/parameters.csv
sed -i 's;goolf-1.7.20;foss-2015b;g' Public_RNA-seq_QC/parameter_files/parameters.csv
sed -i 's;fastqExtension,.fq.gz;fastqExtension,.fastq.gz;' Public_RNA-seq_QC/parameter_files/parameters.csv
sed -i 's;1.102-Java-1.7.0_80;1.119-Java-1.7.0_80;' Public_RNA-seq_QC/parameter_files/parameters.csv
sed -i 's;onekgGenomeFasta,${resDir}/${genomeBuild}/indices/human_g1k_v${human_g1k_vers}.fasta;onekgGenomeFasta,${resDir}/ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/GRCh38.p5.genome.fa;' Public_RNA-seq_QC/parameter_files/parameters.csv
sed -i 's;genesRefFlat,${resDir}/picard-tools/Ensembl${ensemblVersion}/${genomeLatSpecies}.${genomeGrchBuild}.${ensemblVersion}.refflat;genesRefFlat,${resDir}/ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/gencode.v24.chr_patch_hapl_scaff.annotation.refflat;' Public_RNA-seq_QC/parameter_files/parameters.csv
sed -i 's;rRnaIntervalList,${resDir}//picard-tools/Ensembl${ensemblVersion}/${genomeLatSpecies}.${genomeGrchBuild}.${ensemblVersion}.rrna.interval_list;rRnaIntervalList,${resDir}/ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/gencode.v24.chr_patch_hapl_scaff.annotation.rRNA.interval_list;' Public_RNA-seq_QC/parameter_files/parameters.csv
bash Public_RNA-seq_QC/parameter_files/convert.sh Public_RNA-seq_QC/parameter_files/parameters.csv Public_RNA-seq_QC/parameter_files/parameters.converted.csv

sed -i 's;group,umcg-wijmenga;group,umcg-biogen;' Public_RNA-seq_quantification/parameter_files/parameters.csv
sed -i 's;tmp03;tmp04;' Public_RNA-seq_quantification/parameter_files/parameters.csv
sed -i 's;projects/umcg-ndeklein/${project};biogen/input/$1/pipelines/results/;' Public_RNA-seq_quantification/parameter_files/parameters.csv
bash Public_RNA-seq_QC/parameter_files/convert.sh Public_RNA-seq_quantification/parameter_files/parameters.csv Public_RNA-seq_quantification/parameter_files/parameters.converted.csv
echo "DONE"


echo "Make samplesheets"
rm Public_RNA-seq_QC/samplesheet1.csv
mkdir -p Public_RNA-seq_QC/samplesheets/
python make_samplesheet.py sample_annotation/UMG-Target_ALS_RNA_Clinical_Data_06122018.txt Public_RNA-seq_QC/samplesheets/samplesheet_$1_RNA.
python make_samplesheet.py sample_annotation/RNA_Metadata_TALS_2_5July2018.txt Public_RNA-seq_QC/samplesheets/samplesheet_$1_RNA.samplesheet5july2018_
echo "DONE"

echo "Changing prepare script, then running"
sed -i 's;/path/to/molgenis-pipelines/compute5/Public_RNA-seq_QC/;;' Public_RNA-seq_QC/prepare_Public_RNA-seq_QC.sh
sed -i 's;/groups/umcg-wijmenga/tmp04/umcg-ndeklein/molgenis-pipelines/compute5/Public_RNA-seq_QC/;;' Public_RNA-seq_QC/prepare_Public_RNA-seq_QC.sh
sed -i 's;parameters.converted.csv;/groups/umcg-biogen/tmp04/biogen/input/$1/Public_RNA-seq_QC/parameter_files/parameters.converted.csv;' Public_RNA-seq_QC/prepare_Public_RNA-seq_QC.sh
sed -i 's;workflow.csv;/groups/umcg-biogen/tmp04/biogen/input/$1/Public_RNA-seq_QC/workflows/workflowSTAR.csv;' Public_RNA-seq_QC/prepare_Public_RNA-seq_QC.sh
sed -i 's;-rundir /groups/umcg-wijmenga/tmp04/umcg-ndeklein/rundirs/QC/;-rundir /groups/umcg-biogen/tmp04/biogen/input/$1/pipelines/alignment/alignmentDir;' Public_RNA-seq_QC/prepare_Public_RNA-seq_QC.sh
sed -i 's;/groups/umcg-wijmenga/tmp04/umcg-ndeklein/samplesheets/;/groups/umcg-biogen/tmp04/biogen/input/$1/Public_RNA-seq_QC/samplesheets/;' Public_RNA-seq_QC/prepare_Public_RNA-seq_QC.sh

mkdir Public_RNA-seq_QC/prepare_scripts/
for samplesheet in `ls Public_RNA-seq_QC/samplesheets/samplesheet*txt`;
do
    name=$(echo ${samplesheet} | awk -F"." '{print $2}')
    echo "making prepare script using $samplesheet: Public_RNA-seq_QC/prepare_scripts/prepare_Public_RNA-seq_QC.$name.sh"
    rsync -vP Public_RNA-seq_QC/prepare_Public_RNA-seq_QC.sh Public_RNA-seq_QC/prepare_scripts/prepare_Public_RNA-seq_QC.$name.sh
    sed -i "s;samplesheet.csv;samplesheet_$1_RNA.$name.txt;" Public_RNA-seq_QC/prepare_scripts/prepare_Public_RNA-seq_QC.$name.sh
    sed -i "s;alignmentDir;$name;" Public_RNA-seq_QC/prepare_scripts/prepare_Public_RNA-seq_QC.$name.sh
done

cd Public_RNA-seq_QC/prepare_scripts;
for f in *sh;
do
    bash $f;
done
cd -

rm -f Public_RNA-seq_quantification/samplesheet.combined.txt
echo "combining samplesheets for sailfish"
for samplesheet in `ls Public_RNA-seq_QC/samplesheets/samplesheet*txt`;
do
    echo "$samplesheet"
    cat $samplesheet >> Public_RNA-seq_quantification/samplesheet.combined.txt
done
sed -i 's;-rundir results/;-rundir /groups/umcg-biogen/tmp04/biogen/input/$1/pipelines/quantification/;' Public_RNA-seq_quantification/prepare_quantification.sh
sed -i 's;samplesheet.csv;samplesheet.combined.txt;' Public_RNA-seq_quantification/prepare_quantification.sh
cd Public_RNA-seq_quantification/
bash prepare_quantification.sh

cd -
echo "DONE"
