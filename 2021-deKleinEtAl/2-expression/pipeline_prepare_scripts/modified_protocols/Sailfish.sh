#MOLGENIS nodes=1 ppn=1 mem=8gb walltime=23:59:00

### variables to help adding to database (have to use weave)
#string internalId
#string sampleName
#string project
#string uniqueID
###
#string stage
#string checkStage
#string reads1FqGz
#string reads2FqGz
#string sailfishDir
#string sailfishVersion
#string uniqueID
#string sailfishIndex
#string libType
#string numBootstraps
#string flags
#string cramFileDir
#string samtoolsVersion
#string bedtoolsVersion
#string picardVersion
#string iolibVersion
#string onekgGenomeFasta

#Load modules
${stage} Sailfish/${sailfishVersion}
${stage} picard/${picardVersion}
${stage} io_lib/${iolibVersion}
${stage} SAMtools/${samtoolsVersion}

#check modules
${checkStage}

echo "## "$(date)" Start $0"
echo "ID (internalId-project-sampleName): ${internalId}-${project}-${sampleName}"

outDir=${sailfishDir}/${uniqueID}/
mkdir -p $outDir

scramble \
    -I cram \
    -O bam \
    -m \
    -r ${onekgGenomeFasta} \
    ${cramFileDir}${uniqueID}.cram \
    $TMPDIR/${uniqueID}.bam \
    -t 1

returnCode=$?
echo "returncode: $returnCode";
if [ $returnCode -eq 0 ]
then
  echo "convert CRAM to BAM successful"
else
  echo "ERROR: couldn't convert to CRAM to BAM"
  exit 1;
fi


echo "Starting BAM to FASTQ conversion: sort BAM file";
samtools sort \
    -n \
    -o $TMPDIR/${uniqueID}.sorted.bam \
    $TMPDIR/${uniqueID}.bam
rm $TMPDIR/${uniqueID}.bam

returnCode=$?
echo "returncode: $returnCode";
if [ $returnCode -eq 0 ]
then
  echo "sorting BAM successful"
else
  echo "ERROR: couldn't sort BAM"
  exit 1;
fi

# get the filenames from the full path
fq1NameGz=$uniqueID.R1.fq.gz
fq1Name=$uniqueID.R1.fq
fq2NameGz=$uniqueID.R2.fq.gz
fq2Name=$uniqueID.R2.fq


echo "Starting BAM to FASTQ conversion: convert sorted BAM file to fastq"
if [ ${#reads2FqGz} -eq 0 ];
then
  samtools fastq \
      -0 $TMPDIR/$fq1Name \
      $TMPDIR/${uniqueID}.sorted.bam
else
  samtools fastq \
      -1 $TMPDIR/$fq1Name \
      -2 $TMPDIR/$fq2Name \
      $TMPDIR/${uniqueID}.sorted.bam
fi

returnCode=$?
echo "returncode: $returnCode";
if [ $returnCode -eq 0 ]
then
  echo "convert BAM to fastq successful"
else
  echo "ERROR: couldn't convert BAM to fastq"
  exit 1;
fi



if [ ${#reads2FqGz} -eq 0 ]; then
  echo "Single end Sailfish of ${reads1FqGz}"
  # don't know how to determine libType from the fastq files, so defined in parameter file..
  # TODO: add a check if the libtype is compatible with the quant option
  if sailfish quant \
        -i ${sailfishIndex} \
        -l ${libType} \
        -r $TMPDIR/$fq1Name \
        -o $outDir} \
        --numBootstraps ${numBootstraps} \
        ${flags}
  then
    echo "returncode: $?";
  else
    echo "returncode: $?";
    echo "fail";
    exit 1;
  fi
else 
  if sailfish quant \
        -i ${sailfishIndex} \
        -l ${libType} \
        -1 $TMPDIR/$fq1Name -2 $TMPDIR/$fq2Name \
        -o $outDir \
        --numBootstraps ${numBootstraps} \
        ${flags}
  then
    echo "returncode: $?";
  else
    echo "returncode: $?";
    echo "fail";
    exit 1;
  fi
fi

echo "## "$(date)" ##  $0 Done "
