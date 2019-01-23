#MOLGENIS nodes=1 ppn=1 mem=16gb walltime=23:59:00

### variables to help adding to database (have to use weave)
#string sampleName
#string project
###
#list unfilteredBam
#list mergeBamFilesBam
#list sortedBam
#list markDuplicatesBam
#list splitAndTrimBam
#list bqsrBam
echo "## "$(date)" Start $0"

for bam in \${unfilteredBam[@]};do echo "rm ${bam}"; rm ${bam};done

for bam in \${mergeBamFilesBam[@]};do echo "rm ${bam}"; rm ${bam};done

for bam in \${sortedBam[@]};do echo "rm ${bam}"; rm ${bam};done

for bam in \${markDuplicatesBam[@]};do echo "rm ${bam}"; rm ${bam};done

for bam in \${splitAndTrimBam[@]};do echo "rm ${bam}"; rm ${bam};done

for bam in \${bqsrBam[@]};do echo "rm ${bam}"; rm ${bam};done

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
