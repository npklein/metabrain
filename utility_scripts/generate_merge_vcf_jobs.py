
import math
import glob  
import argparse

parser = argparse.ArgumentParser(description='Create slurm job to merge gVCFs')
parser.add_argument('list_of_vcfs', help='.txt file with list of path of VCFs to merge (only for chr15, this will be replaced for other chr numbers)')
parser.add_argument('jobs_dir', help='directory to write jobs to')
parser.add_argument('ref_genome', help='path to reference genome')
parser.add_argument('outdir', help='directory to write merged gVCFs to')
parser.add_argument('project', help='Name of the project (will be prepended to the output file)')
parser.add_argument('chr', help='Chromosome that is being merged')

args = parser.parse_args()

# Below is the job template that will be used
template = """#!/bin/bash
#SBATCH --job-name=MergeGvcfs_chrREPLACECHROMOSOME
#SBATCH --output=MergeGvcfs_chrREPLACECHROMOSOME.out
#SBATCH --error=MergeGvcfs_chrREPLACECHROMOSOME.err
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
    --output REPLACEOUTDIR/REPLACEPROJECT.chrREPLACECHROMOSOME.g.vcf.gz \\
    -L chrREPLACECHROMOSOME \\
    REPLACEINPUT

if [ ! -f REPLACEOUTDIR/REPLACEPROJECT.chrREPLACECHROMOSOME.g.vcf.gz ]; then
    echo "REPLACEOUTDIR/REPLACEPROJECT.chrREPLACECHROMOSOME.g.vcf.gz"
    exit 1
fi

if [ "$?" -eq 0 ];
then
  cd REPLACEOUTDIR/
  md5sum REPLACEPROJECT.chrREPLACECHROMOSOME.g.vcf.gz > REPLACEPROJECT.chrREPLACECHROMOSOME.g.vcf.gz.md5
  cd -
  echo "succes moving files";
else
 echo "returncode: $?";
 echo "fail";
fi


touch REPLACEFINISHED

echo "## "$(date)" ##  $0 Done "


"""

to_merge = ''
# Read in the VCFs to merge
with open(args.list_of_vcfs) as input_file:
    for line in input_file
        if len(line.strip()) == 0:
            continue
        to_merge += ' --variant '+line.strip()

outfile = args.jobs_dir+'MergeGvcfs_chr'+chr+'.sh'
with open(outfile,'w') as out:
    new_template = template.replace('REPLACECHROMOSOME',args.chr)
    new_template = new_template.replace('REPLACEOUTDIR',args.outdir)
    new_template = new_template.replace('REPLACEINPUT', to_merge)
    new_template = new_template.replace('REPLACEPROJECT',args.project)
    new_template = new_template.replace('REPLACEREFGENOME',args.ref_genome)
    new_template = new_template.replace('REPLACEFINISHED','MergeGvcfs_chr'+args.chr+'.sh.finished')
    out.write(new_template)
    print(outfile)

