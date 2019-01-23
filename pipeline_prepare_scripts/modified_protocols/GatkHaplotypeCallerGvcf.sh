#MOLGENIS walltime=23:59:00 mem=14gb ppn=1

### variables to help adding to database (have to use weave)
#string sampleName
#string project
###
#string stage
#string checkStage

#string WORKDIR
#string projectDir
#string gatkVersion
#string dbsnpVcf
#string dbsnpVcfIdx
#string onekgGenomeFasta
#string genomeBuild
#string resDir
#string referenceFastaName
#string bqsrBam
#string chromosome
#string haplotyperDir
#string toolDir
#string haplotyperTargetsPrefix
#string haplotyperTargetsPostfix

echo "## "$(date)" Start $0"

#Load gatk module
${stage} GATK/${gatkVersion}
${checkStage}

mkdir -p ${haplotyperDir}

#do variant calling for all 25 chromosomes seperate
#for this purpose haplotyperGvcf variable is split in ${haplotyperDir}${sampleName}.chr$CHR.g.vcf

$EBROOTGATK/gatk HaplotypeCaller \
   -R ${onekgGenomeFasta} \
   --dbsnp ${dbsnpVcf} \
   -I ${bqsrBam} \
   --dont-use-soft-clipped-bases \
   --standard-min-confidence-threshold-for-calling 10.0 \
   -O ${haplotyperDir}${sampleName}.${chromosome}.g.vcf.gz \
   -L ${chromosome} \
   --emit-ref-confidence GVCF \
   --TMP_DIR=$TMPDIR

returnCode=$?
echo "returncode: $returnCode";
if [ $returnCode -eq 0 ]
then
    echo "Gatk succesful"
else
    echo "ERROR: GATK not successful"
    exit 1;
fi

echo "returncode: $?";
#haplotyperGvcf is split into seperate variables now

cd ${haplotyperDir}
sum $(basename ${haplotyperDir}${sampleName}.${chromosome}.g.vcf.gz)> $(basename ${haplotyperDir}${sampleName}.${chromosome}.g.vcf.gz).md5sum
cd -



echo "## "$(date)" ##  $0 Done "
