#MOLGENIS nodes=1 ppn=1 mem=16gb walltime=5:59:00

### variables to help adding to database (have to use weave)
#string sampleName
#string project
###
#string stage
#string checkStage
#string samtoolsVersion
#string gatkVersion
#string markDuplicatesBam
#string markDuplicatesBai
#string onekgGenomeFasta
#string chromosome

#string splitAndTrimBam
#string splitAndTrimBai
#string splitAndTrimDir
#string toolDir

echo "## "$(date)" Start $0"



${stage} SAMtools/${samtoolsVersion}
${stage} GATK/${gatkVersion}
${checkStage}

mkdir -p ${splitAndTrimDir}

$EBROOTGATK/gatk SplitNCigarReads \
 -R ${onekgGenomeFasta} \
 -I ${markDuplicatesBam} \
 -O ${splitAndTrimBam} \
 --TMP_DIR=$TMPDIR \
 -L ${chromosome}

returnCode=$?
echo "returncode: $returnCode";
if [ $returnCode -eq 0 ]
then
  echo "SplitAndTrim succesful"
else
  echo "ERROR: Could not SplitAndTrim"
  exit 1;
fi


echo "returncode: $?";

echo "## "$(date)" ##  $0 Done "
