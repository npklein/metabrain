#MOLGENIS nodes=1 ppn=1 mem=8gb walltime=23:59:00

### variables to help adding to database (have to use weave)
#string internalId
#string sampleName
#string project
#string uniqueID
###
#string stage
#string checkStage
#string brainDir
#string enaSamplesheet
#string reads1FqGz

#Load modules
${stage} Python/${pythonVersion}

#check modules
${checkStage}

echo "## "$(date)" Start $0"
echo "ID (internalId-project-sampleName): ${internalId}-${project}-${sampleName}"


python $brainDir/ENA/download_ENA_samples.py ${enaSamplesheet} \
                                             $(dirname ${reads1FqGz}) \
                                             --sample ${sampleName}

returnCode=$?
echo "returncode: $returnCode";
if [ $returnCode -eq 0 ]
then
  echo "Downloading fastq file successful"
else
  echo "ERROR: Could not download fastq file"
  exit 1;
fi



echo "## "$(date)" ##  $0 Done "
