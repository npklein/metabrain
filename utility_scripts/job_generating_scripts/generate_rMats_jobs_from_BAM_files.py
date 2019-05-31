# making rmats jobs directly from bam files, don't need to convert cram to bam like original
import argparse
import glob
import os

parser = argparse.ArgumentParser(description='Generate rMATs files.')
parser.add_argument('output_directory', help='Outputdir to write results to')
parser.add_argument('jobs_directory', help='Directory to write jobs to')
parser.add_argument('bam_base_directory', help='Directory with bamFiles')
parser.add_argument('libblas_location', help='Location of directory with libblas library')
parser.add_argument('gtf', help='Location of gtf file')
parser.add_argument('--sample_split', help='Character to split file name on and take [0] as sample name', default='.cram')

args = parser.parse_args()

bam_files = glob.glob(args.bam_base_directory+'/**/*bam', recursive=True)
print('found ',len(bam_files),'bam files')

outdir = args.output_directory
job_base_dir = args.jobs_directory

if not os.path.exists(outdir):
    os.makedirs(outdir)
if not os.path.exists(job_base_dir):
    os.makedirs(job_base_dir)

def make_jobs(template):
    for bam in bam_files:
        if 'TCX' in bam:
            batch = 'MayoTCX'
        else:
            batch = bam.split('/')[-2]
        jobs_dir = job_base_dir + '/'+batch+'/'
        if not os.path.exists(jobs_dir):
            os.mkdir(jobs_dir)
        
        if batch == 'MSBB':
            sample = bam.split('/')[-1].split('.accepted_hits')[0].split('_resequenced')[0]
        elif batch == 'MayoCBE':
            sample = bam.split('/')[-1].split('.')[0]
        elif batch == 'MayoTCX' or batch == 'ROSMAP':
            sample = bam.split('/')[-1].split('.bam')[0].split('Aligned')[0]
        else:
            print(bam)
            raise RuntimeError('Unknown batch: '+batch)

        results_dir = outdir+'/'+batch
        if not os.path.exists(results_dir):
            os.mkdir(results_dir)

        new_template = template.replace('REPLACENAME', sample)
        new_template = new_template.replace('REPLACEOUT', results_dir+'/'+sample)
        new_template = new_template.replace('REPLACEBAMCOPY', bam.replace('.bam','.copy.bam'))
        new_template = new_template.replace('REPLACEBAM', bam.replace('bam','bam'))
        new_template = new_template.replace('REPLACELIBBLAS',args.libblas_location)
        new_template = new_template.replace('REPLACEGTF',args.gtf)
        with open(jobs_dir+'/'+sample+'.sh','w') as out:
            out.write(new_template)



template = '''#!/bin/bash
#SBATCH --job-name=rMats_REPLACENAME
#SBATCH --output=REPLACENAME.out
#SBATCH --error=REPLACENAME.err
#SBATCH --time=05:59:00
#SBATCH --cpus-per-task 1
#SBATCH --mem 8gb
#SBATCH --nodes 1
module load rMATS/4.0.2-foss-2015b-Python-2.7.11
module load SAMtools

rsync -vP REPLACEBAM $TMPDIR/$(basename REPLACEBAMCOPY)

echo "REPLACEBAM" > $TMPDIR/b1.txt
echo "$TMPDIR/$(basename REPLACEBAMCOPY)" > $TMPDIR/b2.txt

# test if it is single end or paired end.
# 1. get the header out of the BAM file (samtools view -H). This contains the use command that was used to make BAM file
# 2. get the line out of the header that contains the user command (grep)
# 3. split the line on "user command line" (awk) and print everything that is after
# 4. in this string, count the number of occurences of fastq.gz and fq.gz
FASTQINPUTFILE=$(samtools view -H REPLACEBAM | grep "user command line" | awk -F"user command line:" '{ print $2}' | grep -o "\.fastq.gz\|\.fq.gz\|\.fastq" | wc -l)

t=
if [ "$FASTQINPUTFILE" -eq 1 ];
then
    t="single"
    echo "BAM file is single-end, will only extract 1 fastq file"
elif [ "$FASTQINPUTFILE" -eq 2 ];
then
    t="paired"
    echo "BAM file is paired-end, will extract 2 fastq files"
else
    echo "ERROR: PAIRED was $PAIRED, should have been 1 or 2";
    exit 1;
fi

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:REPLACELIBBLAS/
python $EBROOTRMATS/rMATS-turbo-Linux-UCS4/rmats.py \\
    --b1 $TMPDIR/b1.txt \\
    --b2 $TMPDIR/b2.txt \\
    --gtf REPLACEGTF \\
    --od REPLACEOUT \\
    -t $t \\
    --readLength 100 \\
    --cstat 0.0001 \\
    --libType fr-unstranded

if [ $? -eq 0 ];
then
    echo "succes!"
    echo "returncode: $?"
    rm $TMPDIR/$(basename REPLACEBAMCOPY)
else
    echo "FAIL!"
    echo "returncode: $?"
    exit 1;
fi

'''


make_jobs(template)
