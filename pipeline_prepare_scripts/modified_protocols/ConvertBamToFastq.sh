#MOLGENIS nodes=1 ppn=4 mem=12gb walltime=05:59:00

### variables to help adding to database (have to use weave)
#string internalId
#string sampleName
#string project
###
#string stage
#string checkStage
#string picardVersion
#string iolibVersion
#string reads1FqGz
#string reads2FqGz
#string samtoolsVersion
#string alignedBam

# Get input file

echo "IMPORTANT!! This works only for paired-end samples"

#Load modules
${stage} picard/${picardVersion}
${stage} io_lib/${iolibVersion}
${stage} SAMtools/${samtoolsVersion}

#check modules
${checkStage}

echo "## "$(date)" Start $0"
echo "ID (internalId-project-sampleName): ${internalId}-${project}-${sampleName}"


echo "Starting BAM to FASTQ conversion: sort BAM file";
samtools sort \
    -@ 4 \
    -n \
    -o $TMPDIR/$(basename $alignedBam) \
    ${alignedBam}

returnCode=$?
echo "returncode: $returnCode";
if [ $returnCode -eq 0 ]
then
  echo "sorting BAM successful"
else
  echo "ERROR: couldn't sort BAM"
  exit 1;
fi

samtools fastq \
  -@ 4 \
  -1 $TMPDIR/$(basename ${reads1FqGz%.gz}) \
  -2 $TMPDIR/$(basename ${reads2FqGz%.gz}) \
  $TMPDIR/$(basename $alignedBam)

echo "gzipping:"
echo "gzip $TMPDIR/$(basename ${reads1FqGz%.gz})"
echo "gzip $TMPDIR/$(basename ${reads2FqGz%.gz})"
gzip $TMPDIR/$(basename ${reads1FqGz%.gz})
gzip $TMPDIR/$(basename ${reads2FqGz%.gz})



rsync -vP $TMPDIR/$(basename ${reads1FqGz}) ${reads1FqGz}
rsync -vP $TMPDIR/$(basename ${reads2FqGz}) ${reads2FqGz}
returnCode=$?

echo "returncode: $returnCode";
if [ $returnCode -eq 0 ]
then
  echo "convert BAM to fastq successful"
else
  echo "ERROR: couldn't convert BAM to fastq"
  exit 1;
fi


if [ $returnCode -eq 0 ]
then
  echo "merge successful, remove $alignedBam"
  rm $alignedBam
else
  echo "ERROR: couldn't merge fastq files"
  exit 1;
fi


echo "final returncode: $?";
echo "succes removing files";


echo "## "$(date)" ##  $0 Done "






