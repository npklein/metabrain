import argparse
import glob
import os

parser = argparse.ArgumentParser(description='Generate rMATs files.')
parser.add_argument('cram_base_directory', help='Directory with cramFiles')
parser.add_argument('jobs_directory', help='Directory to write jobs to')
parser.add_argument('output_directory', help='Outputdir to write results to')
parser.add_argument('ref_gtf', help='Reference gtf file location')
parser.add_argument('feature_type', help='Feature type to count (e.g. exon or transcript)')
parser.add_argument('--extra_options', help='Extra options to give to featurecounts (e.g. -O --fraction',
                     nargs='+')


args = parser.parse_args()
extra_options = []
if args.extra_options:
    for x in args.extra_options:
        if len(x) == 1:
            x = '-'+x
        else:
            x = '--'+x
        extra_options.append(x)
    extra_options = ' '.join(extra_options)
else:
    extra_options = ''

cram_files = glob.glob(args.cram_base_directory+'/**/*cram', recursive=True)+glob.glob(args.cram_base_directory+'/**/*bam', recursive=True)
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
        new_template = new_template.replace('REPLACEBAM', cram.replace('cram','bam'))
        new_template = new_template.replace('REPLACEGTF',args.ref_gtf)
        new_template = new_template.replace('REPLACEFEATURETYPE',args.feature_type)
        new_template = new_template.replace('REPLACEEXTRAOPTIONS',extra_options)
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

featureCounts \
    -f \ # read summarization will be performed at featurelevel  (eg.   exon  level)
    -C \ # chimeric fragments (those fragments that have their  two  ends  aligned  to  different  chromosomes)  will  NOT be counted
    -s 0 \ # un-stranded
    -p \ #  fragments (or templates) will be counted insteadof reads.  This option is only applicable for paired-end reads.
    -t REPLACEFEATURETYPE \ #  feature type.  Only rows which have the matchedfeature type in the provided GTF annotation file will be in-cluded for read counting
    -g gene_id \ # attribute type used to group features (eg. exons) into meta-features (eg. genes) when GTF annotation is pro-vided. ‘geneid’ by default. This attribute type is usually thegene identifier.
    -O \ # reads (or fragments if-pis specified) will be al-lowed to be assigned to more than one matched meta-feature
    -−fraction \ # each overlapping feature will receive a count of 1/y, where y is the total number of features overlapping with the read
    -a REPLACEGTF \
    -o REPLACEOUT \
    $TMPDIR/$(basename REPLACEBAM) \
    REPLACEEXTRAOPTIONS

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
