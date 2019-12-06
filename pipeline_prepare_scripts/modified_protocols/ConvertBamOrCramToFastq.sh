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
#string alignedBamOrCram

# Get input file

mkdir -p $(dirname $reads1FqGz)

#Load modules
${stage} picard/${picardVersion}
${stage} io_lib/${iolibVersion}
${stage} SAMtools/${samtoolsVersion}

#check modules
${checkStage}

echo "## "$(date)" Start $0"
echo "ID (internalId-project-sampleName): ${internalId}-${project}-${sampleName}"


echo "Starting BAM/CRAM to FASTQ conversion: sort BAM file";

samtools sort \
  -@ 4 \
  -n \
  -o $TMPDIR/$(basename ${alignedBamOrCram}) \
  ${alignedBamOrCram}

returnCode=$?
echo "returncode: $returnCode";
if [ $returnCode -eq 0 ]
then
  echo "sorting BAM/CRAM successful"
else
  echo "ERROR: couldn't sort BAM/CRAM"
  exit 1;
fi

if [[ -n "$reads2FqGz" ]];
then
  echo "paired-end, writing 2 fastq files"
  samtools fastq \
    -@ 4 \
    -1 $TMPDIR/$(basename ${reads1FqGz}) \
    -2 $TMPDIR/$(basename ${reads2FqGz}) \
    $TMPDIR/$(basename $alignedBamOrCram)
    rsync -vP $TMPDIR/$(basename ${reads1FqGz}) ${reads1FqGz}
    rsync -vP $TMPDIR/$(basename ${reads2FqGz}) ${reads2FqGz}
else
  echo "single-end, writing 1 fastq file"
  samtools fastq \
    -@ 4 \
    -0 $TMPDIR/$(basename ${reads1FqGz}) \
    $TMPDIR/$(basename $alignedBamOrCram)
    rsync -vP $TMPDIR/$(basename ${reads1FqGz}) ${reads1FqGz}
fi

returnCode=$?

echo "returncode: $returnCode";
if [ $returnCode -eq 0 ]
then
  echo "convert BAM/CRAM to fastq successful"
else
  echo "ERROR: couldn't convert BAM/CRAM to fastq"
  exit 1;
fi



if [[ ! -f "${reads1FqGz}" ]];
then
    echo "ERROR: ${reads1FqGz} does not exist"
    exit 1;
fi

if [[ ${#reads2FqGz} -eq 0 && -f "${reads2FqGz}" ]];
then
    echo "ERROR: ${reads2FqGz} does not exist"
fi


echo "final returncode: $?";
echo "succes removing files";


echo "## "$(date)" ##  $0 Done "






