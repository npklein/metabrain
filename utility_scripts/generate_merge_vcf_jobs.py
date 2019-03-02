
import math
import glob  
import argparse

parser = argparse.ArgumentParser(description='Create slurm jobs to merge gVCFs together in batches of 200')
parser.add_argument('list_of_vcfs', help='.txt file with list of path of VCFs to merge (only for chr15, this will be replaced for other chr numbers)')
parser.add_argument('jobs_dir', help='directory to write jobs to')
parser.add_argument('ref_genome', help='path to reference genome')
parser.add_argument('outdir', help='directory to write merged gVCFs to')

args = parser.parse_args()

# Below is the job template that will be used
# For each chromosome samples are divided up in 200 and for each batch of 200, one job is made
template = """#!/bin/bash
#SBATCH --job-name=MergeGvcfs_batchREPLACEBATCH_chrREPLACECHROMOSOME
#SBATCH --output=MergeGvcfs_batchREPLACEBATCH_chrREPLACECHROMOSOME.out
#SBATCH --error=MergeGvcfs_batchREPLACEBATCH_chrREPLACECHROMOSOME.err
#SBATCH --time=23:59:59
#SBATCH --cpus-per-task 1
#SBATCH --mem 18gb
#SBATCH --nodes 1
#SBATCH --export=NONE
#SBATCH --get-user-env=L

ENVIRONMENT_DIR="."
set -e
set -u

echo "## "$(date)" Start $0"

#Load gatk module
module load GATK/4.0.8.1-foss-2015b-Python-3.6.3
module list

mkdir -p REPLACEOUTDIR


$EBROOTGATK/gatk CombineGVCFs \\
    --reference REPLACEREFGENOME \\
    --output REPLACEOUTDIR/REPLACEPROJECT.batchREPLACEBATCH_chrREPLACECHROMOSOME.g.vcf.gz \\
    -L chrREPLACECHROMOSOME \\
    REPLACEINPUT

then
 echo "returncode: $?"; 

if [ ! -f REPLACEOUTDIR/REPLACEPROJECT.batchREPLACEBATCH_chrREPLACECHROMOSOME.g.vcf.gz ]; then
    echo "REPLACEOUTDIR/REPLACEPROJECT.batchREPLACEBATCH_chrREPLACECHROMOSOME.g.vcf.gz"
    exit 1
fi
cd REPLACEOUTDIR/
md5sum REPLACEPROJECT.batchREPLACEBATCH_chrREPLACECHROMOSOME.g.vcf.gz > REPLACEPROJECT.batchREPLACEBATCH_chrREPLACECHROMOSOME.g.vcf.gz.md5
 cd -
 echo "succes moving files";
else
 echo "returncode: $?";
 echo "fail";
fi


touch REPLACEFINISHED

echo "## "$(date)" ##  $0 Done "


"""

# Read in the VCFs to merge
with open(args.list_of_vcfs) as input:
    samples = 0
    lines = input.read().split('\n')

# Count number of samples (skipping empty lines)
for line in lines:
    if len(line.strip()) == 0:
        continue
    samples += 1

# Make batch per chromosome
batches = math.ceil(samples/200)
print('writing to',batches,'batches')
batch = 0
input = {'1':{},'2':{},'3':{},'4':{},'5':{},'6':{},'7':{},
         '8':{},'9':{},'10':{},'11':{},'12':{},'13':{},'14':{},
         '15':{},'16':{},'17':{},'18':{},'19':{},'20':{},
         '21':{},'22':{},'X':{},'Y':{},'M':{}}


lines = sorted(lines)
for line in lines:
    if len(line.strip()) == 0:
        continue
    for chr in range(1, 26, 1):
        if chr == 23:
            chr = 'X'
        elif chr == 24:
            chr = 'Y'
        elif chr == 25:
            chr = 'M'
        if str(batch) in input[str(chr)]:
            input[str(chr)][str(batch)] += ' --variant '+line.strip().replace('chr15','chr'+str(chr))
        else:
            input[str(chr)][str(batch)] = ' --variant '+line.strip().replace('chr15','chr'+str(chr))            
    batch += 1
    if batch >= batches:
        batch = 0

for chr in input:
    for batch in input[chr]:
        outfile = args.jobs_dir+'MergeGvcfs_batch'+batch+'_chr'+chr+'.sh'
        with open(outfile,'w') as out:
            new_template = template.replace('REPLACEBATCH',batch)
            new_template = new_template.replace('REPLACECHROMOSOME',chr)
            new_template = new_template.replace('REPLACEOUTDIR',args.outdir)
            new_template = new_template.replace('REPLACEINPUT', input[chr][batch])
            new_template = new_template.replace('REPLACEPROJECT','BIOS')
            new_template = new_template.replace('REPLACEREFGENOME',args.ref_genome)
            new_template = new_template.replace('REPLACEFINISHED','MergeGvcfs_batch'+batch+'_chr'+chr+'.sh.finished')
            out.write(new_template)
        print(outfile)

