import argparse
import glob
import os

parser = argparse.ArgumentParser(description='Generate rMATs files.')
parser.add_argument('output_directory', help='Outputdir to write results to')
parser.add_argument('jobs_directory', help='Directory to write jobs to')
parser.add_argument('cram_base_directory', help='Directory with cramFiles')
parser.add_argument('libblas_location', help='Location of directory with libblas library')

args = parser.parse_args()

cram_files = glob.glob('/groups/umcg-biogen/tmp04/input/rawdata/*/pipelines/results/cramFiles/*cram', recursive=True)
print('found ',len(cram_files),'cram files')

outdir = '/groups/umcg-biogen/tmp04/umcg-ndeklein/rMats/rMats_results/'
job_base_dir = 'rMats_jobs/'

if not os.path.exists(outdir):
    raise RuntimeError('output directory does not exist')

def make_jobs(template):
    x = 0
    batch_number = 0
    jobs_dir = job_base_dir + 'batch'+str(batch_number)+'/'
    if not os.path.exists(jobs_dir):
        os.mkdir(jobs_dir)
    for cram in cram_files:
        x += 1
        if x % 25 == 0:
            batch_number += 1
            jobs_dir = job_base_dir + 'batch'+str(batch_number)+'/'
            if not os.path.exists(jobs_dir):
                os.mkdir(jobs_dir)

        sample = cram.split('/')[-1].split('.cram')[0]
        new_template = template.replace('REPLACENAME', sample)
        new_template = new_template.replace('REPLACEOUT', outdir+'/'+sample)
        new_template = new_template.replace('REPLACECRAM', cram)
        new_template = new_template.replace('REPLACEBAMCOPY', cram.replace('.cram','.copy.bam'))
        new_template = new_template.replace('REPLACEBAM', cram.replace('cram','bam'))
        new_template = new_template.replace('REPLACELIBBLAS',args.libblas_location)
        with open(jobs_dir+'/'+sample+'.sh','w') as out:
            out.write(new_template)



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
    rm $TMPDIR/$(basename REPLACEBAM)
    rm $TMPDIR/$(basename REPLACEBAMCOPY)
else
    echo "FAIL!"
    echo "returncode: $?"
    exit 1;
fi

'''


make_jobs(template)
