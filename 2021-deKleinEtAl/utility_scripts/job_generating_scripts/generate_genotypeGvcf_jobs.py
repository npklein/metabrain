import os
import math
import glob  
import argparse

parser = argparse.ArgumentParser(description='Create slurm jobs to call genotype on a genomicDB')
parser.add_argument('genomicDB_dir', help='Directory that contains the merged gVCFs')
parser.add_argument('jobs_dir', help='directory to write jobs to')
parser.add_argument('ref_genome', help='path to reference genome')
parser.add_argument('outdir', help='directory to write merged gVCFs to')
parser.add_argument('dbsnp', help='dbSNP file')
parser.add_argument('project', help='Name of the project (will be prepended to the file name)')
parser.add_argument('intervallist_dir', help='Directory with GATK curated interval lists, split into chrName of the project (file name format: hg38.even.handcurated.20k.chr<..>.intervals)')

args = parser.parse_args()

template = """#!/bin/bash
#SBATCH --job-name=GenotypeGvcf_REPLACEINTERVALNAME
#SBATCH --output=GenotypeGvcf_REPLACEINTERVALNAME.out
#SBATCH --error=GenotypeGvcf_REPLACEINTERVALNAME.err
#SBATCH --time=23:59:59
#SBATCH --cpus-per-task 1
#SBATCH --mem 20gb
#SBATCH --nodes 1
#SBATCH --export=NONE
#SBATCH --get-user-env=L


ENVIRONMENT_DIR="."
set -e
set -u

echo "## "$(date)" Start $0"



#Load gatk module
module load GATK/4.1.0.0-foss-2015b-Python-3.6.3
module list

mkdir -p REPLACEOUTDIR


$EBROOTGATK/gatk --java-options "-Xmx20G" GenotypeGVCFs \\
    --reference REPLACEREFGENOME \\
    --dbsnp REPLACEDBSNP \\
    --output $TMPDIR/REPLACEOUTPUT \\
    --variant gendb://REPLACEINPREFIX.REPLACECHROMOSOME \\
    --stand-call-conf 10.0 \\
    -L REPLACEINTERVAL \\
    -G StandardAnnotation 

rsync -rvP $TMPDIR/REPLACEOUTPUT* REPLACEOUTDIR/
if [ "$?" -eq 0 ];
then
 echo "returncode: $?"; 
 
 cd REPLACEOUTDIR
 md5sum REPLACEOUTPUT > REPLACEOUTPUT.md5
 cd -
 echo "succes moving files";
else
 echo "returncode: $?";
 echo "fail";
fi
touch REPLACESHSCRIPT.finished
echo "## "$(date)" ##  $0 Done "I

"""

chromosomes = ['1','2','3','4','5','6','7',
        '8','9','10','11','12','13','14',
         '15','16','17','18','19','20',
         '21','22','X','Y']  
# There are no curated intervals for chrM, if interested have to call on whole chr
#,'M']

for chr in chromosomes:
    chr = 'chr'+chr
    print(chr)
    if not os.path.exists(args.outdir+'/'+chr):
        os.makedirs(args.outdir+'/'+chr)
    if not os.path.exists(args.jobs_dir+'/'+chr):
        os.makedirs(args.jobs_dir+'/'+chr)
    
    with open(args.intervallist_dir+'/hg38.even.handcurated.20k.'+chr+'.intervals') as input_file:
        for interval in input_file:
            interval = interval.strip()
            interval_name = interval.replace(':','_').replace('-','_')
            outfile = args.jobs_dir+'/'+chr+'/GenotypeGvcf_'+interval_name+'.sh'
            with open(outfile,'w') as out:
                new_template = template.replace('REPLACECHROMOSOME',chr)
                new_template = new_template.replace('REPLACEOUTDIR',args.outdir+'/'+chr+'/')
                new_template = new_template.replace('REPLACEOUTPUT',args.project+'_'+interval_name+'.gg.vcf.gz')
                new_template = new_template.replace('REPLACEINPREFIX', args.genomicDB_dir+'/'+args.project)
                new_template = new_template.replace('REPLACEREFGENOME',args.ref_genome)
                new_template = new_template.replace('REPLACEDBSNP', args.dbsnp)
                new_template = new_template.replace('REPLACEINTERVALNAME', interval_name)
                new_template = new_template.replace('REPLACEINTERVAL', interval)
                new_template = new_template.replace('REPLACESHSCRIPT', args.jobs_dir+'/'+chr+'/GenotypeGvcf_'+interval_name+'.sh')
                out.write(new_template)

