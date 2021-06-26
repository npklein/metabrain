import argparse
import glob
import os

parser = argparse.ArgumentParser(description='Generate htseq job files.')
parser.add_argument('jobs_directory', help='Directory for jobs')
parser.add_argument('cram_base_directory', help='Directory with cramFiles')
parser.add_argument('ref_gtf', help='Reference gtf file location')
parser.add_argument('--sample_split', help='Character to split file name on and take [0] as sample name', default='.cram')


args = parser.parse_args()

cram_files = glob.glob(args.cram_base_directory+'/**/*cram', recursive=True)
print('found ',len(cram_files),'cram files')

job_base_dir = args.jobs_directory

if not os.path.exists(job_base_dir):
    os.mkdir(job_base_dir)

def make_jobs(template):
    batch_number = 0
    x = 0
    studies = set([])
    for cram in cram_files:
        assert(cram.endswith('.cram'))
        if x % 25 == 0:
            batch_number += 1
        # /scratch/umcg-ndeklein/tmp03/ENA/pipelines/no_patch_chromosomes/results/PRJNA231202_batch0/cramFiles/SAMN02441931_SRR1047827.cram
        
        jobs_dir = job_base_dir+'/batch'+str(batch_number)
        if not os.path.exists(jobs_dir):
            os.mkdir(jobs_dir)

        sample = cram.split('/')[-1].split(args.sample_split)[0]
        new_template = template.replace('REPLACENAME', sample)
        new_template = new_template.replace('REPLACEOUT', cram.split('/cramFiles/')[0]+'/star/'+sample+'.htSeqCount.txt')
        new_template = new_template.replace('REPLACECRAM', cram)
        new_template = new_template.replace('REPLACEBAM', cram.replace('.cram','.bam'))
        new_template = new_template.replace('REPLACEGTF',args.ref_gtf)
        with open(jobs_dir+'/'+sample+'.sh','w') as out:
            out.write(new_template)
        x += 1
    print('Added: '+','.join(studies))


template = '''#!/bin/bash
#SBATCH --job-name=HtSeq_REPLACENAME
#SBATCH --output=REPLACENAME.out
#SBATCH --error=REPLACENAME.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task 1
#SBATCH --mem 8gb
#SBATCH --nodes 1
module load HTSeq/0.6.1p1-foss-2015b-Python-2.7.10
module load SAMtools

echo $TMPDIR
echo "converting cram to bam"
samtools view -hb REPLACECRAM > $TMPDIR/$(basename REPLACEBAM)

htseq-count -f bam \\
     $TMPDIR/$(basename REPLACEBAM) \\
     REPLACEGTF \\
    >> REPLACEOUT

if [ $? -eq 0 ];
then
    echo "succes!"
    echo "returncode: $?"
    rm $TMPDIR/$(basename REPLACEBAM)
else
    echo "FAIL!"
    echo "returncode: $?"
    exit 1;
fi

'''


make_jobs(template)
