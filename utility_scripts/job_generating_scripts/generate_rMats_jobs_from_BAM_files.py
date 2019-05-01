# making rmats jobs directly from bam files, don't need to convert cram to bam like original
import argparse
import glob
import os

parser = argparse.ArgumentParser(description='Generate rMATs files.')
parser.add_argument('output_directory', help='Outputdir to write results to')
parser.add_argument('jobs_directory', help='Directory to write jobs to')
parser.add_argument('bam_base_directory', help='Directory with bamFiles')
parser.add_argument('libblas_location', help='Location of directory with libblas library')

args = parser.parse_args()

bam_files = glob.glob(args.bam_base_directory+'/**/*bam', recursive=True)
print('found ',len(bam_files),'bam files')

outdir = args.output_directory
job_base_dir = args.jobs_directory

if not os.path.exists(outdir):
    os.mkdir(outdir)
if not os.path.exists(job_base_dir):
    os.mkdir(job_base_dir)

def make_jobs(template):
    for bam in bam_files:
        batch = bam.split('/')[-2]
        jobs_dir = job_base_dir + '/'+batch+'/'
        if not os.path.exists(jobs_dir):
            os.mkdir(jobs_dir)

        sample = bam.split('/')[-1].split('.accepted_hits')[0].split('_resequenced')[0]
        new_template = template.replace('REPLACENAME', sample)
        new_template = new_template.replace('REPLACEOUT', outdir+'/'+sample)
        new_template = new_template.replace('REPLACEBAMCOPY', bam.replace('.bam','.copy.bam'))
        new_template = new_template.replace('REPLACEBAM', bam.replace('bam','bam'))
        new_template = new_template.replace('REPLACELIBBLAS',args.libblas_location)
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

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:REPLACELIBBLAS/
python $EBROOTRMATS/rMATS-turbo-Linux-UCS4/rmats.py \\
    --b1 $TMPDIR/b1.txt \\
    --b2 $TMPDIR/b2.txt \\
    --gtf /apps/data/ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/gencode.v24.chr_patch_hapl_scaff.annotation.gtf \\
    --od REPLACEOUT \\
    -t paired \\
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
