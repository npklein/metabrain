#MOLGENIS nodes=1 ppn=2 mem=14gb walltime=23:59:00

### variables to help adding to database (have to use weave)
#string sampleName
#string project
###
#string stage
#string checkStage
#string samtoolsVersion
#string gatkVersion
#string onekgGenomeFasta
#string goldStandardVcf
#string goldStandardVcfIdx
#string oneKgPhase1IndelsVcf
#string oneKgPhase1IndelsVcfIdx
#string dbsnpVcf
#string dbsnpVcfIdx
#string indelRealignmentBam
#string indelRealignmentBai
#string bqsrDir
#string bqsrBam
#string bqsrBai
#string bqsrBeforeGrp
#string toolDir

#pseudo from gatk forum (link: http://gatkforums.broadinstitute.org/discussion/3891/best-practices-for-variant-calling-on-rnaseq):
#java -jar GenomeAnalysisTK.jar -T SplitNCigarReads -R ref.fasta -I dedupped.bam -o split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

echo "## "$(date)" Start $0"



${stage} GATK/${gatkVersion}
${checkStage}


mkdir -p ${bqsrDir}

#do bqsr for covariable determination then do print reads for valid bqsrbams
#check the bqsr part and add known variants

$EBROOTGATK/gatk BaseRecalibrator\
 -R ${onekgGenomeFasta} \
 -I ${indelRealignmentBam} \
 -o ${bqsrBeforeGrp} \
 -knownSites ${dbsnpVcf} \
 -knownSites ${goldStandardVcf}\
 -knownSites ${oneKgPhase1IndelsVcf}\
 -nct 2 \
 --TMP_DIR=$TMPDIR

returnCode=$?
echo "returncode: $returnCode";
if [ $returnCode -eq 0 ]
then
  echo "BaseCalibrator succesful"
else
  echo "ERROR: Could not calibrate bases"
  exit 1;
fi


$EBROOTGATK/gatk \
 -T PrintReads \
 -R ${onekgGenomeFasta} \
 -I ${indelRealignmentBam} \
 -o ${bqsrBam} \
 -BQSR ${bqsrBeforeGrp} \
 -nct 2 \
 --TMP_DIR=$TMPDIR

returnCode=$?
echo "returncode: $returnCode";
if [ $returnCode -eq 0 ]
then
  echo "PrintReads succesful"
else
  echo "ERROR: Could not PrintReads"
  exit 1;
fi



echo "returncode: $?";

cd ${bqsrDir}
md5sum $(basename ${bqsrBam})> $(basename ${bqsrBam}).md5sum
md5sum $(basename ${bqsrBam%bam}bai)> $(basename ${bqsrBam%bam}bai).md5sum
cd -

echo "## "$(date)" ##  $0 Done "
