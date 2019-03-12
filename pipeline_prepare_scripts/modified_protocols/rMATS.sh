#MOLGENIS nodes=1 ppn=1 mem=8gb walltime=06:00:00

### variables to help adding to database (have to use weave)
#string internalId
#string sampleName
#string project
###
#string sortedBamDir
#string stage
#string checkStage
#string picardVersion
#string iolibVersion
#string sortedBamDir
#string cramFileDir
#string onekgGenomeFasta
#string uniqueID
#string rmatsVersion

module load rMATS/${rmatsVersion}
module load SAMtools/${samtoolsVersion}
echo "converting cram to bam"
samtools view -hb ${cramFileDir}${uniqueID}.cram > $TMPDIR/$(basename ${cramFileDir}${uniqueID}.bam)
samtools view -hb ${cramFileDir}${uniqueID}.cram > $TMPDIR/$(basename ${cramFileDir}${uniqueID}.bam)

echo "DONE"
echo "$(basename /scratch/umcg-ndeklein/tmp03/ENA/pipelines/results/PRJNA212047_batch2/bamFiles/SAMN02251119_SRR934792.bam)" > b1.txt
echo "$(basename /scratch/umcg-ndeklein/tmp03/ENA/pipelines/results/PRJNA212047_batch2/bamFiles/SAMN02251124_SRR934801.bam)" > b2.txt

#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scratch/umcg-ndeklein/tmp03/ENA/rMats/rMATS-turbo-Linux-UCS4/usr/lib64/
echo "start rMATS"
python /scratch/umcg-ndeklein/tmp03/ENA/rMats/rMATS-turbo-Linux-UCS4/rmats.py \
    --b1 b1.txt \
    --b2 b2.txt \
    --gtf /data/umcg-ndeklein/apps/data/ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/gencode.v24.chr_patch_hapl_scaff.annotation.gtf \
    --od /scratch/umcg-ndeklein/tmp03/ENA/rMats/rMats_results/SAMN02251119_SRR934792___SAMN02251124_SRR934801 \
    -t paired \
    --readLength 100 \
    --cstat 0.0001 \
    --libType fr-unstranded

for f in /scratch/umcg-ndeklein/tmp03/ENA/rMats/rMats_results/SAMN02251119_SRR934792___SAMN02251124_SRR934801/*JCEC.txt;
do
    if [[ ! $(wc -l <$f) -ge 2 ]]
    then
        echo "ERROR: $f has max 1 line, rMATS did not run correctly"
        exit 1;
    fi
done

if [ $? -eq 0 ];
then
    echo "succes!"
    echo "returncode: $?"
    rm $TMPDIR/$(basename /scratch/umcg-ndeklein/tmp03/ENA/pipelines/results/PRJNA212047_batch2/bamFiles/SAMN02251119_SRR934792.bam)
    rm $TMPDIR/$(basename /scratch/umcg-ndeklein/tmp03/ENA/pipelines/results/PRJNA212047_batch2/bamFiles/SAMN02251124_SRR934801.bam)
else
    echo "FAIL!"
    echo "returncode: $?"
    exit 1;
fi

