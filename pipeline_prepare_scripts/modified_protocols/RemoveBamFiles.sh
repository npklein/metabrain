#MOLGENIS nodes=1 ppn=1 mem=16gb walltime=23:59:00

### variables to help adding to database (have to use weave)
#string project
###
#string unfilteredBamDir
#string mergeBamFilesDir
#string sortedBamDir
#string markDuplicatesDir
#string splitAndTrimDir
#string bqsrDir
echo "## "$(date)" Start $0"

if [ -z "$unfilteredBamDir" ];
then
    >&2 echo "ERROR: $unfilteredBamDir is empty!"
    exit 1;
fi
rm ${unfilteredBamDir}/*bam

if [ -z "$mergeBamFilesDir" ];
then
    >&2 echo "ERROR: $mergeBamFilesDir is empty!"
    exit 1;
fi
rm ${mergeBamFilesDir}/*bam

if [ -z "$sortedBamDir" ];
then
    >&2 echo "ERROR: $sortedBamDir is empty!"
    exit 1;
fi
rm ${sortedBamDir}/*bam

if [ -z "$markDuplicatesDir" ];
then
    >&2 echo "ERROR: $markDuplicatesDir is empty!"
    exit 1;
fi
rm ${markDuplicatesDir}/*bam

if [ -z "$splitAndTrimDir" ];
then
    >&2 echo "ERROR: $splitAndTrimDir is empty!"
    exit 1;
fi
rm ${splitAndTrimDir}/*bam

if [ -z "$bqsrDir" ];
then
    >&2 echo "ERROR: $bqsrDir is empty!"
    exit 1;
fi
rm ${bqsrDir}/*bam

returnCode=$?
echo "returncode: $returnCode";
if [ $returnCode -eq 0 ]
then
  echo "Removing BAMs succesful"
else
  echo "ERROR: Could not remove BAMs"
  exit 1;
fi



echo "## "$(date)" ##  $0 Done "
