#MOLGENIS nodes=1 ppn=1 mem=14gb walltime=1:59:00

### variables to help adding to database (have to use weave)
#string sampleName
#string project
###
#string stage
#string checkStage
#string RVersion
#string gatkVersion
#string onekgGenomeFasta
#string goldStandardVcf
#string goldStandardVcfIdx
#string oneKgPhase1IndelsVcf
#string oneKgPhase1IndelsVcfIdx
#string dbsnpVcf
#string dbsnpVcfIdx
#string bqsrDir
#string bqsrBam
#string bqsrBai
#string analyseCovarsDir
#string bqsrBeforeGrp
#string bqsrAfterGrp
#string analyseCovariatesPdf
#string toolDir
#string analyseCovariatesIntermediateCsv

echo "## "$(date)" Start $0"



${stage} R/${RVersion}
${stage} GATK/${gatkVersion}
${checkStage}

mkdir -p ${analyseCovarsDir}

#do bqsr for covariable determination then do print reads for valid bqsrbams
#check the bqsr part and add known variants

$EBROOTGATK/gatk BaseRecalibrator\
 -R ${onekgGenomeFasta} \
 -I ${bqsrBam} \
 -O ${bqsrAfterGrp} \
 --known-sites ${dbsnpVcf} \
 --known-sites ${goldStandardVcf} \
 --known-sites ${oneKgPhase1IndelsVcf} \
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

$EBROOTGATK/gatk  AnalyzeCovariates \
 --ignore-last-modification-times \
 --before-report-file ${bqsrBeforeGrp} \
 --after-report-file ${bqsrAfterGrp} \
 --verbosity DEBUG \
 --intermediate-csv-file ${analyseCovariatesIntermediateCsv} \
 --plots-report-file ${analyseCovariatesPdf} \
 --TMP_DIR=$TMPDIR

returnCode=$?
echo "returncode: $returnCode";
if [ $returnCode -eq 0 ]
then
  echo "AnalyzeCovariates succesful"
else
  echo "ERROR: Could not analyze covariates"
  exit 1;
fi

echo "returncode: $?";

echo "## "$(date)" ##  $0 Done "
