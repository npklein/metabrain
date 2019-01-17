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
#list bqsrBam
#string haplotyperDir
#string toolDir
#string haplotyperTargetsPrefix
#string haplotyperTargetsPostfix

echo "## "$(date)" Start $0"

#Load gatk module
${stage} GATK/${gatkVersion}
${checkStage}

#sort unique and print like 'INPUT=file1.bam INPUT=file2.bam '
bams=($(printf '%s\n' "${bqsrBam[@]}" | sort -u ))

inputs=$(printf ' -I %s ' $(printf '%s\n' ${bams[@]}))

mkdir -p ${haplotyperDir}

#do variant calling for all 25 chromosomes seperate
#for this purpose haplotyperGvcf variable is split in ${haplotyperDir}${sampleName}.chr$CHR.g.vcf

for CHR in {1..25}
do
   echo "CHR $CHR"
   $EBROOTGATK/gatk HaplotypeCaller \
       -R ${onekgGenomeFasta} \
       --dbsnp ${dbsnpVcf} \
       $inputs \
       -dontUseSoftClippedBases \
       -stand_call_conf 10.0 \
       -stand_emit_conf 20.0 \
       -O ${haplotyperDir}${sampleName}.chr$CHR.g.vcf.gz \
       -variant_index_type LINEAR \
       -variant_index_parameter 128000 \
       -L ${haplotyperTargetsPrefix}$CHR${haplotyperTargetsPostfix} \
       --emitRefConfidence GVCF \
       --TMP_DIR=$TMPDIR

       returnCode=$?
       echo "returncode: $returnCode";
       if [ $returnCode -eq 0 ]
       then
          echo "Gatk succesful for chr $CHR"
       else
          echo "ERROR: GATK not successful for chr $CHR"
          exit 1;
       fi

    echo "returncode: $?";
    #haplotyperGvcf is split into seperate variables now

    cd ${haplotyperDir}
	md5sum $(basename ${haplotyperDir}${sampleName}.chr$CHR.g.vcf.gz)> $(basename ${haplotyperDir}${sampleName}.chr$CHR.g.vcf.gz).md5sum
    cd -

done


echo "## "$(date)" ##  $0 Done "
