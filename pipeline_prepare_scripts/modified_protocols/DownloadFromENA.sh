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

#Load modules
${stage} Python/${pythonVersion}

#check modules
${checkStage}

echo "## "$(date)" Start $0"
echo "ID (internalId-project-sampleName): ${internalId}-${project}-${sampleName}"


python $brainDir/ENA/download_ENA_samples.py

returnCode=$?
echo "returncode: $returnCode";
if [ $returnCode -eq 0 ]
then
  echo "convert CRAM to BAM successful"
else
  echo "ERROR: couldn't convert to CRAM to BAM"
  exit 1;
fi



echo "## "$(date)" ##  $0 Done "
