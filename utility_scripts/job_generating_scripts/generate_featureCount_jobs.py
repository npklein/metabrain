import argparse
import glob
import os

parser = argparse.ArgumentParser(description='Generate featureCounts files.')
parser.add_argument('cram_base_directory', help='Directory with cramFiles')
parser.add_argument('jobs_directory', help='Directory to write jobs to')
parser.add_argument('output_directory', help='Outputdir to write results to')
parser.add_argument('ref_gtf', help='Reference gtf file location')
parser.add_argument('ref_meta_gtf', help='Reference gtf file location')
#parser.add_argument('--extra_options', help='Extra options to give to featurecounts (e.g. -O --fraction',
#                     nargs='+')
#parser.add_argument('feature_type', help='Feature type to count (e.g. exon or transcript)')


args = parser.parse_args()
#extra_options = []
#if args.extra_options:
#    for x in args.extra_options:
#        if len(x) == 1:
#            x = '-'+x
#        else:
#            x = '--'+x
#        extra_options.append(x)
#    extra_options = ' '.join(extra_options)
#else:
#    extra_options = ''

cram_files = glob.glob(args.cram_base_directory+'/**/*.cram', recursive=True)+glob.glob(args.cram_base_directory+'/**/*.bam', recursive=True)
print('found ',len(cram_files),'cram and bam files')

outdir = args.output_directory
job_base_dir = args.jobs_directory

if not os.path.exists(outdir):
    os.makedirs(outdir)
if not os.path.exists(job_base_dir):
    os.makedirs(job_base_dir)

def make_jobs(template):
    prev_study = None
    for cram in cram_files:
        if not cram.endswith('.cram') and not cram.endswith('.bam'):
            continue
        if 'AMP_AD' in cram:
            study = cram.split('/')[-2]
        elif 'CMC_HBCC' in cram:
            study = 'CMC_HBCC'
        else:
            study = cram.split('/pipelines/')[0].split('/')[-2]
        if study == 'BPD':
            continue
        if study == '':
            print(cram)
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
        elif study == 'ucl-upload-biogen':
            sample = cram.split('/')[-1].split('.')[0]
            study = 'Braineac'    
        elif study == 'psychEncode' or study == 'CMC' or study == 'CMC_HBCC':
            sample = cram.split('/')[-1].split(".cram")[0].replace("individualID.","").replace("specimenID.","")
            if '.Aligned' in sample:
                sample = sample.split('.Aligned')[0]
        else:
            print(study)
            print(cram)
            raise RuntimeError('Unknown study: '+study)
        new_template = template.replace('REPLACENAME', sample)
        new_template = new_template.replace('REPLACEOUT', outdir+'/'+study+'/')
        new_template = new_template.replace('REPLACECRAM', cram)
        new_template = new_template.replace('REPLACEBAM', cram.replace('cram','bam'))
        new_template = new_template.replace('REPLACEGTF',args.ref_gtf)
        new_template = new_template.replace('REPLACEMETAEXONGTF',args.ref_meta_gtf)
#        new_template = new_template.replace('REPLACEFEATURETYPE',args.feature_type)
#        new_template = new_template.replace('REPLACEEXTRAOPTIONS',extra_options)
        with open(jobs_dir+'/'+sample+'.sh','w') as out:
            out.write(new_template)
        prev_study = study


template = '''#!/bin/bash
#SBATCH --job-name=featureCounts_REPLACENAME
#SBATCH --output=REPLACENAME.out
#SBATCH --error=REPLACENAME.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task 1
#SBATCH --mem 8gb
#SBATCH --nodes 1

set -e
set -u

module load Subread/1.6.4-foss-2015b
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

echo "Using $INPUTBAM as input file"

#   -f          read summarization will be performed at featurelevel  (eg.   exon  level)
#   -C          chimeric fragments (those fragments that have their  two  ends  aligned  to  different  chromosomes)  will  NOT be counted
#   -s 0        un-stranded
#   -p          fragments (or templates) will be counted insteadof reads.  This option is only applicable for paired-end reads.
#   -t exon     feature type.  Only rows which have the matchedfeature type in the provided GTF annotation file will be in-cluded for read counting
#   -g gene_id  attribute type used to group features (eg. exons) into meta-features (eg. genes) when GTF annotation is pro-vided. ‘geneid’ by default. This attribute type is usually thegene identifier.
#   -O          reads (or fragments if-pis specified) will be al-lowed to be assigned to more than one matched meta-feature
#   --fraction  each overlapping feature will receive a count of 1/y, where y is the total number of features overlapping with the read

echo "Start exon.countAll"
mkdir -p REPLACEOUT/exon.countAll
TMPOUT=$TMPDIR/REPLACENAME.exon.countAll.txt
featureCounts -f -C -s 0 -p -t exon -g gene_id -O \\
    -a REPLACEGTF \\
    -o $TMPOUT \\
    $INPUTBAM
gzip $TMPOUT
rsync -vP $TMPOUT* REPLACEOUT/exon.countAll/

echo "Start exon.countFraction"
mkdir -p REPLACEOUT/exon.countFraction
TMPOUT=$TMPDIR/REPLACENAME.exon.countFraction.txt
featureCounts -f -C -s 0 -p -t exon -g gene_id -O --fraction \\
    -a REPLACEGTF \\
    -o $TMPOUT \\
    $INPUTBAM
gzip $TMPOUT
rsync -vP $TMPOUT* REPLACEOUT/exon.countFraction/

echo "Start metaExon.countAll"
mkdir -p REPLACEOUT/metaExon.countAll
TMPOUT=$TMPDIR/REPLACENAME.metaExon.countAll.txt
featureCounts -f -C -s 0 -p -t exonic_part -g gene_id -O \\
    -a REPLACEMETAEXONGTF \\
    -o $TMPOUT \\
    $INPUTBAM
gzip $TMPOUT
rsync -vP $TMPOUT* REPLACEOUT/metaExon.countAll/

echo "Start metaExon.countFraction"
mkdir -p REPLACEOUT/metaExon.countFraction
TMPOUT=$TMPDIR/REPLACENAME.metaExon.countFraction.txt
featureCounts -f -C -s 0 -p -t exonic_part -g gene_id -O --fraction \\
    -a REPLACEMETAEXONGTF \\
    -o $TMPOUT \\
    $INPUTBAM
gzip $TMPOUT
rsync -vP $TMPOUT* REPLACEOUT/metaExon.countFraction/

echo "Start transcript.countAll"
mkdir -p REPLACEOUT/transcript.countAll/
TMPOUT=$TMPDIR/REPLACENAME.transcript.countAll.txt
featureCounts -f -C -s 0 -p -t transcript -g gene_id -O \\
    -a REPLACEGTF \\
    -o $TMPOUT \\
    $INPUTBAM
gzip $TMPOUT
rsync -vP $TMPOUT* REPLACEOUT/transcript.countAll//

echo "Start transcript.countFraction"
mkdir -p REPLACEOUT/transcript.countFraction
TMPOUT=$TMPDIR/REPLACENAME.transcript.countFraction.txt
featureCounts -f -C -s 0 -p -t transcript -g gene_id -O --fraction \\
    -a REPLACEGTF \\
    -o $TMPOUT \\
    $INPUTBAM

gzip $TMPOUT
rsync -vP $TMPOUT* REPLACEOUT/transcript.countFraction/

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
