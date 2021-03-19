#MOLGENIS walltime=24:00:00 nodes=1 cores=1 mem=6gb

#Parameter mapping
#string onekgGenomeFasta
#string cramFileDir
#string annotationGtf
#string htseqTxtOutput
#string samtoolsVersion
#string htseqVersion
#string stranded
#string htseqDir
#string mode
#string sortType
#string featureType
#string internalId
#string sampleName
#string project
#string uniqueID
#string pysamVersion
#string unfilteredTwoPassBamDir
#string projectDir
#string iolibVersion
#string uniqueID

echo "## "$(date)" Start $0"
echo "ID (internalId-project-sampleName): ${internalId}-${project}-${sampleName}"
#Echo parameter values
bam="${cramFileDir}${uniqueID}.cram"
annotationGtf="${annotationGtf}"
htseqTxtOutput="${htseqTxtOutput}"

echo -e "bam=${cramFileDir}${uniqueID}.cram\nannotationGtf=${annotationGtf}\nhtseqTxtOutput=${htseqTxtOutput}"

module load SAMtools/${samtoolsVersion}
module load Pysam/${pysamVersion}
module load HTSeq/${htseqVersion}
module load io_lib/${iolibVersion}
module list

echo "Assuming that the bam file is position sorted, if htseq fails check if your input bam is sorted"
mkdir -p ${htseqDir}

echo -e "\nQuantifying expression"

mkdir -p ${projectDir}/htseqTwoPass/
#samtools sort -o ${unfilteredTwoPassBamDir}/${uniqueID}.sorted.bam ${unfilteredTwoPassBamDir}/${uniqueID}.bam
scramble \
    -I cram \
    -O bam \
    -m \
    -r ${onekgGenomeFasta} \
    ${cramFileDir}${uniqueID}.cram \
    $TMPDIR/${uniqueID}.bam \
    -t 1


python /home/umcg-ndeklein/brain_eQTL/pipeline_prepare_scripts/external_scripts//dexseq_count.py -p yes -s no -r pos -f bam \
        /apps/data/ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/exonicGFF/gencode.v24.chr_patch_hapl_scaff.annotation.DEXSeq.gff \
         $TMPDIR/${uniqueID}.bam ${htseqTxtOutput%txt}exonic.txt
#if htseq-count \
#        -m ${mode} \
#        -f bam \
#        -t ${featureType} \
#        --stranded ${stranded} \
#        ${unfilteredTwoPassBamDir}/${uniqueID}.sorted.bam \
#        -r pos \
#        ${annotationGtf} >  ${htseqTxtOutput}___tmp___ ;
#then
#        if [[ $(wc -l <${htseqTxtOutput}___tmp___) -ge 2 ]]
#        then
#            echo "Gene count succesfull"
#            echo "returncode: $?"
#            mv ${htseqTxtOutput}___tmp___ ${projectDir}/htseqTwoPass/$(basename ${htseqTxtOutput})
#        else
#            echo "output not written correctly";
 #           echo "returncode: 1"
#            exit 1;
#        fi
#else
#        echo "Genecount failed"
#        echo 'returncode: 1'
#        exit 1
#fi

echo "## "$(date)" ##  $0 Done "
