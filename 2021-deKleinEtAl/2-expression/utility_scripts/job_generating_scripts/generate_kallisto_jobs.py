import sys
import argparse
import glob
import os

parser = argparse.ArgumentParser(description='Generate featureCounts files.')
parser.add_argument('cram_base_directory', help='Directory with cramFiles')
parser.add_argument('jobs_directory', help='Directory to write jobs to')
parser.add_argument('output_directory', help='Outputdir to write results to')
parser.add_argument('kallisto_index', help='Kallisto index file location')
parser.add_argument('reference', help='Reference file that was used to make cram file')
parser.add_argument('-k','--kallisto_version', help='Version of kallisto to use',default='0.43.1-foss-2015b')
parser.add_argument('-q','--qos', help='QOS to use',default='regular')
parser.add_argument('-s','--string_to_ignore', help='If string in filename, ignore the file. e.g. -p AMP-AD ignores all files with AMP-AD in file path')


args = parser.parse_args()

cram_files = glob.glob(args.cram_base_directory+'/**/*.cram', recursive=True)#+glob.glob(args.cram_base_directory+'/**/*.bam', recursive=True)
print('exluding dirs with no_patch_chromosomes in it')
cram_files = [x for x in cram_files if 'no_patch_chromosomes' in x]

if args.string_to_ignore:
    print('Removing files with string '+args.string_to_ignore)
    n_cramfiles = len(cram_files)
    cram_files = [x for x in cram_files if args.string_to_ignore not in x]
    print(n_cramfiles - len(cram_files),'files removed because of string: '+args.string_to_ignore)

print('found ',len(cram_files),'cram and bam files')
sys.stdout.flush()

unique_project = set(['/'.join(x.split('/')[:-1]) for x in cram_files])
print('Unique directories where cramfiles where taken from:')
print(unique_project)
print('----')
outdir = args.output_directory
job_base_dir = args.jobs_directory

if not os.path.exists(outdir):
    os.makedirs(outdir)
if not os.path.exists(job_base_dir):
    os.makedirs(job_base_dir)

def make_jobs(cram_files, template):
    prev_study = None
    for cram in cram_files:
        if not (cram.endswith('.cram') and not cram.endswith('.bam')) or '/BPD/' in cram:
            continue
        study = None
        if 'AMP_AD' in cram:
            if 'MayoTCX' in cram:
                study='MayoTCX'
            elif 'MSBB' in cram:
                study='MSBB'
            elif 'MayoCBE' in cram:
                study='MayoCBE'
            elif 'ROSMAP' in cram:
                study='ROSMAP'
        elif 'CMC_HBCC' in cram:
            study = 'CMC_HBCC'
        elif 'BipSeq' in cram:
            study = 'BipSeq'
        elif 'UCLA_ASD' in cram:
            study = 'UCLA_ASD'
        elif 'BrainGVEx' in cram:
            study = 'BrainGVEx'
        elif 'GTEx' in cram:
            study = 'GTEx'
        elif 'TargetALS' in cram:
            study = 'TargetALS'
        elif 'MayoTCX' in cram:
            study = 'MayoTCX'
        elif 'MayoCBE' in cram:
            study = 'MayoCBE'
        elif 'ROSMAP' in cram:
            study = 'ROSMAP'
        elif 'ENA' in cram:
            study = 'ENA'
        elif 'CMC' in cram:
            study = 'CMC'
        elif 'NABEC' in cram:
            study = 'NABEC'
        elif 'ucl-upload-biogen' in cram:
            study = 'Braineac'
        elif 'Brainseq' in cram:
            study = 'Brainseq'
        if not study:
            print(cram)
            raise RuntimeError('Study not set')

        jobs_dir = job_base_dir + '/'+study+'/'
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
        elif study == 'Braineac':
            sample = cram.split('/')[-1].split('.')[0]
            study = 'Braineac'    
        elif study == 'BipSeq' or study == 'CMC' or study == 'CMC_HBCC' or study == 'BrainGVEx' or study == 'CMC_HBCC' or study == 'UCLA_ASD':
            sample = cram.split('/')[-1].split(".cram")[0].replace("individualID.","").replace("specimenID.","")
        else:
            print(study)
            sys.stdout.flush()
            print(cram)
            sys.stdout.flush()
            raise RuntimeError('Unknown study: '+study)
        new_template = template.replace('REPLACENAME', sample)
        new_template = new_template.replace('REPLACEOUT', outdir+'/'+study+'/'+sample+'/')
        new_template = new_template.replace('REPLACECRAM', cram)
        new_template = new_template.replace('REPLACEBAM', cram.replace('cram','bam'))
        new_template = new_template.replace('REPLACEINDEX',args.kallisto_index)
        new_template = new_template.replace('REPLACEKALLISTOVERSION',args.kallisto_version)
        new_template = new_template.replace('REPLACEREFERENCE',args.reference)
        new_template = new_template.replace('REPLACEQOS',args.qos)

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
#SBATCH --qos=REPLACEQOS

set -e
set -u

module load kallisto/REPLACEKALLISTOVERSION
module load SAMtools

if [[ "$REPLACECRAM" == *cram ]];
then
    echo "converting cram to bam"
    samtools view -hb \\
                  -T REPLACEREFERENCE \\
                   REPLACECRAM > $TMPDIR/$(basename REPLACEBAM)
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
  echo "sorting BAM finished"
else
  echo "ERROR: couldn't sort BAM"
  exit 1;
fi

TMPFASTQ1=$TMPDIR/$(basename ${INPUTBAM%.bam}.R1.fastq)
TMPFASTQ2=$TMPDIR/$(basename ${INPUTBAM%.bam}.R2.fastq)
TMPFASTQ3=$TMPDIR/$(basename ${INPUTBAM%.bam}.rest.fastq)

# test if it is single end or paired end. 
# 1. get the header out of the BAM file (samtools view -H). This contains the use command that was used to make BAM file
# 2. get the line out of the header that contains the user command (grep)
# 3. split the line on "user command line" (awk) and print everything that is after
# 4. in this string, count the number of occurences of fastq.gz and fq.gz
FASTQINPUTFILE=$(samtools view -H $TMPDIR/$(basename $INPUTBAM) | grep "user command line" | awk -F"--readFilesIn" '{ print $2}' | grep -o "\.fastq.gz\|\.fq.gz\|\.fastq" | wc -l)

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
    echo "FASTQUINPUTFILE was $FASTQINPUTFILE"
    exit 1;
fi


TMPOUT=$TMPDIR/REPLACENAME/
if [ "$PAIRED" -eq 1 ];
then
    samtools fastq \\
        -@ 4 \\
        -1 $TMPFASTQ1 \\
        -2 $TMPFASTQ2 \\
        $TMPDIR/$(basename $INPUTBAM) > $TMPFASTQ3

    kallisto quant \\
        --bias \\
        -t 4 \\
        -i REPLACEINDEX \\
        -o $TMPOUT \\
        $TMPFASTQ1 \\
        $TMPFASTQ2 

elif [ "$PAIRED" -eq 0 ];
then
    samtools fastq \\
        -@ 4 \\
        -1 $TMPFASTQ1 \\
        $TMPDIR/$(basename $INPUTBAM) > $TMPFASTQ3
    
    # see if the file size of TMPFASTQ3 > TMPFATQ1, if so, use that
    filesize1=$(stat -c%s $TMPFASTQ1)
    filesize2=$(stat -c%s $TMPFASTQ3)

    if [[ "$filesize1" -gt "$filesize2" ]]; 
    then
        echo "$TMPFASTQ1 larger than $TMPFASTQ3, use for kallisto"
    else
        echo "$TMPFASTQ3 larger than $TMPFASTQ1, use for kallisto"
        TMPFASTQ1=$TMPFASTQ3
    fi
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


make_jobs(cram_files, template)
