import argparse
import glob
import os

parser = argparse.ArgumentParser(description='Generate rMATs files.')
parser.add_argument('output_directory', help='Outputdir to write results to')
parser.add_argument('jobs_directory', help='Directory to write jobs to')
parser.add_argument('cram_base_directory', help='Directory with cramFiles')
parser.add_argument('libblas_location', help='Location of directory with libblas library')
parser.add_argument('ref_gtf', help='Reference gtf file location')
parser.add_argument('--sample_split', help='Character to split file name on and take [0] as sample name', default='.cram')


args = parser.parse_args()

cram_files = glob.glob(args.cram_base_directory+'/**/*cram', recursive=True)
print('found ',len(cram_files),'cram files')

outdir = args.output_directory
job_base_dir = args.jobs_directory

if not os.path.exists(outdir):
    os.mkdir(outdir)
if not os.path.exists(job_base_dir):
    os.mkdir(job_base_dir)

def make_jobs(template):
    x = 0
    studies = set([])
    for cram in cram_files:
        if not cram.endswith('.cram') and not cram.endswith('.bam') or '/BPD/' in cram:
            continue
        study = None
        if 'AMP_AD' in cram:
            study = cram.split('/')[-2]
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
        elif 'CMC' in cram:
            study = 'CMC'
        elif 'NABEC' in cram:
            study = 'NABEC'
        elif 'ucl-upload-biogen' in cram:
            study = 'Braineac'
        elif 'Brainseq' in cram:
            study = 'Brainseq'
        elif 'ENA' in cram:
            study = 'ENA'
        if not study:
            print(cram)
            raise RuntimeError('Study not set')
        studies.add(study)
        jobs_dir = job_base_dir +study+'/' 
        if not os.path.exists(jobs_dir):
            os.mkdir(jobs_dir)
        x += 1
        results_dir = outdir+'/'+study
        if not os.path.exists(results_dir):
            os.mkdir(results_dir)

        sample = cram.split('/')[-1].split(args.sample_split)[0]
        new_template = template.replace('REPLACENAME', sample)
        new_template = new_template.replace('REPLACEOUT', results_dir+'/'+sample)
        new_template = new_template.replace('REPLACECRAM', cram)
        new_template = new_template.replace('REPLACEBAMCOPY', cram.replace('.cram','.copy.bam'))
        new_template = new_template.replace('REPLACEBAM', cram.replace('cram','bam'))
        new_template = new_template.replace('REPLACELIBBLAS',args.libblas_location)
        new_template = new_template.replace('REPLACEGTF',args.ref_gtf)
        with open(jobs_dir+'/'+sample+'.sh','w') as out:
            out.write(new_template)
    print('Added: '+','.join(studies))


template = '''#!/bin/bash
#SBATCH --job-name=rMats_REPLACENAME
#SBATCH --output=REPLACENAME.out
#SBATCH --error=REPLACENAME.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task 1
#SBATCH --mem 8gb
#SBATCH --nodes 1
module load rMATS/4.0.2-foss-2015b-Python-2.7.11
module load SAMtools

echo $TMPDIR
echo "converting cram to bam"
samtools view -hb REPLACECRAM > $TMPDIR/$(basename REPLACEBAM)
rsync -vP $TMPDIR/$(basename REPLACEBAM) $TMPDIR/$(basename REPLACEBAMCOPY)

echo "$TMPDIR/$(basename REPLACEBAM)" > $TMPDIR/b1.txt
echo "$TMPDIR/$(basename REPLACEBAMCOPY)" > $TMPDIR/b2.txt

# test if it is single end or paired end.
# 1. get the header out of the BAM file (samtools view -H). This contains the use command that was used to make BAM file
# 2. get the line out of the header that contains the user command (grep)
# 3. split the line on "user command line" (awk) and print everything that is after
# 4. in this string, count the number of occurences of fastq.gz and fq.gz
FASTQINPUTFILE=$(samtools view -H $REPLACECRAM | grep "user command line" | awk -F"user command line:" '{ print $2}' | grep -o ".fastq.gz\|.fq.gz" | wc -l)

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
    -t ${t} \\
    --readLength 100 \\
    --cstat 0.0001 \\
    --libType fr-unstranded

if [ $? -eq 0 ];
then
    echo "succes!"
    echo "returncode: $?"
    rm $TMPDIR/$(basename REPLACEBAM)
    rm $TMPDIR/$(basename REPLACEBAMCOPY)
else
    echo "FAIL!"
    echo "returncode: $?"
    exit 1;
fi

'''


make_jobs(template)
