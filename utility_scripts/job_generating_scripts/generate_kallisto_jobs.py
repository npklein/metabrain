import argparse
import glob
import os

parser = argparse.ArgumentParser(description='Generate featureCounts files.')
parser.add_argument('cram_base_directory', help='Directory with cramFiles')
parser.add_argument('jobs_directory', help='Directory to write jobs to')
parser.add_argument('output_directory', help='Outputdir to write results to')
parser.add_argument('kallisto_index', help='Kallisto index file location')


args = parser.parse_args()

cram_files = glob.glob(args.cram_base_directory+'/**/*.cram', recursive=True)+glob.glob(args.cram_base_directory+'/**/*.bam', recursive=True)
print('found ',len(cram_files),'cram and bam files')

outdir = args.output_directory
job_base_dir = args.jobs_directory

if not os.path.exists(outdir):
    os.makedirs(outdir)
if not os.path.exists(job_base_dir):
    os.makedirs(job_base_dir)

def make_jobs(template):
    x = 0
    batch_number = 0
    prev_study = None
    for cram in cram_files:
        if not cram.endswith('.cram') and not cram.endswith('.bam'):
            continue
        if 'AMP_AD' in cram:
            study = cram.split('/')[-2]
        elif 'CMC_HBCC' in cram:
            study = 'CMC_HBCC'
        else:
            study = cram.split('/pipelines/')[0].split('/')[-1]
        if study == 'BPD':
            continue
        if study != prev_study:
            batch_number = 0
            x = 0
        if study == '':
            print(cram)
        jobs_dir = job_base_dir + '/'+study+'/batch'+str(batch_number)+'/'
        if not os.path.exists(jobs_dir):
            os.makedirs(jobs_dir)
        x += 1
        if x % 25 == 0:
            batch_number += 1
            jobs_dir = job_base_dir + '/'+study+'/batch'+str(batch_number)+'/'
            if not os.path.exists(jobs_dir):
                os.makedirs(jobs_dir)
        if study == 'MSBB':
            sample = cram.split('/')[-1].split('.accepted_hits')[0].split('_resequenced')[0]
        elif study == 'MayoCBE' or study == 'GTEx' or study == 'TargetALS' or study == 'NABEC' or study == 'Brainseq' or study == 'ENA':
            sample = cram.split('/')[-1].split('.')[0]
        elif study == 'MayoTCX' or study == 'ROSMAP':
            if 'Aligned.out' in cram:
                sample = cram.split('/')[-1].split('Aligned')[0]
            else:
                sample = cram.split('/')[-1].split('.cram')[0]
        elif study == 'ucl-upload-biogen':
            sample = cram.split('/')[-1].split('.')[0]
            study = 'Braineac'    
        elif study == 'psychEncode' or study == 'CMC' or study == 'CMC_HBCC':
            sample = cram.split('/')[-1].split(".cram")[0].replace("individualID.","").replace("specimenID.","")
        else:
            print(study)
            print(cram)
            raise RuntimeError('Unknown study: '+study)
        new_template = template.replace('REPLACENAME', sample)
        new_template = new_template.replace('REPLACEOUT', outdir+'/'+study+'/'+sample+'/')
        new_template = new_template.replace('REPLACECRAM', cram)
        new_template = new_template.replace('REPLACEBAM', cram.replace('cram','bam'))
        new_template = new_template.replace('REPLACEINDEX',args.kallisto_index)
        with open(jobs_dir+'/'+sample+'.sh','w') as out:
            out.write(new_template)
        prev_study = study


template = '''#!/bin/bash
#SBATCH --job-name=kallisto_REPLACENAME
#SBATCH --output=REPLACENAME.out
#SBATCH --error=REPLACENAME.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task 4
#SBATCH --mem 16gb
#SBATCH --nodes 1

set -e
set -u

module load Kallisto/0.43.1-foss-2015b
module load SAMtools

if [[ "$REPLACECRAM" == *cram ]];
then
    echo "converting cram to bam"
    samtools view -hb REPLACECRAM > $TMPDIR/$(basename REPLACEBAM)
    INPUTBAM=$TMPDIR/$(basename REPLACEBAM)
else
    echo "Input file is already in bam format, not conversion needed"
    INPUTBAM=REPLACECRAM
fi

echo "Starting BAM ($INPUTBAM) to FASTQ conversion: sort BAM file";
samtools sort \\
    -@ 4 \\
    -n \\
    -o $TMPDIR/$(basename $INPUTBAM) \\
    ${INPUTBAM}

returnCode=$?
echo "returncode: $returnCode";
if [ $returnCode -eq 0 ]
then
  echo "sorting BAM successful"
else
  echo "ERROR: couldn't sort BAM"
  exit 1;
fi

TMPFASTQ1=$TMPDIR/$(basename ${INPUTBAM%.bam}.R1.fastq)
TMPFASTQ2=$TMPDIR/$(basename ${INPUTBAM%.bam}.R2.fastq)

# test if it is single end or paired end. 
# 1. get the header out of the BAM file (samtools view -H). This contains the use command that was used to make BAM file
# 2. get the line out of the header that contains the user command (grep)
# 3. split the line on "user command line" (awk) and print everything that is after
# 4. in this string, count the number of occurences of fastq.gz and fq.gz
FASTQINPUTFILE=$(samtools view -H $TMPDIR/$(basename $INPUTBAM) | grep "user command line" | awk -F"user command line:" '{ print $2}' | grep -o ".fastq.gz\|.fq.gz" | wc -l)

PAIRED=
if [ "$FASTQINPUTFILE" -eq 1 ];
then
    PAIRED=0
    echo "BAM file is single-end, will only extract 1 fastq file"
elif [ "$FASTQINPUTFILE" -eq 2 ];
then
    PAIRED=1
    echo "BAM file is paired-end, will extract 2 fastq files"
else
    echo "ERROR: PAIRED was $PAIRED, should have been 1 or 2";
    exit 1;
fi


TMPOUT=$TMPDIR/REPLACENAME/
if [ "$PAIRED" -eq 1 ];
then
    samtools fastq \
        -@ 4 \
        -1 $TMPFASTQ1 \\
        -2 $TMPFASTQ2 \\
        $TMPDIR/$(basename $INPUTBAM)

    kallisto quant \\
        --bias \\
        -t 4 \\
        -i REPLACEINDEX \\
        -o $TMPOUT \\
        $TMPFASTQ1 \\
        $TMPFASTQ2 

elif [ "$PAIRED" -eq 0 ];
then
    samtools fastq \
        -@ 4 \
        -1 $TMPFASTQ1 \\
        $TMPDIR/$(basename $INPUTBAM)

    kallisto quant \\
        --bias -t 4 \\
        --single \\
        --fragment-length=200 \\
        --sd=20 \\
        -i REPLACEINDEX \\
        -o $TMPOUT \\
        $TMPFASTQ1
fi

mkdir -p REPLACEOUT
rsync -rvP $TMPOUT REPLACEOUT



if [ $? -eq 0 ];
then
    echo "succes!"
    echo "returncode: $?"
    if [[ "$REPLACECRAM" == *cram ]];
    then
        rm $INPUTBAM
    fi
else
    echo "FAIL!"
    echo "returncode: $?"
    exit 1;
fi

'''


make_jobs(template)
