
import math
import glob  
import argparse

parser = argparse.ArgumentParser(description='Create slurm job to Concatenate gVCFs')
parser.add_argument('list_of_vcfs', help='.txt file with list of path of VCFs to Concatenate (only for chr15, this will be replaced for other chr numbers)')
parser.add_argument('jobs_dir', help='directory to write jobs to')
parser.add_argument('outdir', help='directory to write Concatenated gVCFs to')
parser.add_argument('project', help='Name of the project (will be prepended to the output file)')
parser.add_argument('chr', help='Chromosome that is being Concatenated')

args = parser.parse_args()

# Below is the job template that will be used
template = """#!/bin/bash
#SBATCH --job-name=ConcatenateGvcfs_chrREPLACECHROMOSOME
#SBATCH --output=ConcatenateGvcfs_chrREPLACECHROMOSOME.out
#SBATCH --error=ConcatenateGvcfs_chrREPLACECHROMOSOME.err
#SBATCH --time=5:59:59
#SBATCH --cpus-per-task 8
#SBATCH --mem 20gb
#SBATCH --nodes 1
#SBATCH --export=NONE
#SBATCH --get-user-env=L

ENVIRONMENT_DIR="."
set -e
set -u

echo "## "$(date)" Start $0"

#Load gatk module
module load BCFtools/1.9-foss-2018a
module list

mkdir -p REPLACEOUTDIR


bcftools concat --threads 8 \\
                -o REPLACEOUTDIR/REPLACEPROJECT.chrREPLACECHROMOSOME.gg.vcf.gz \\
                -O z \\
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

to_Concatenate = ''
# Read in the VCFs to Concatenate
with open(args.list_of_vcfs) as input_file:
    for line in input_file:
        if len(line.strip()) == 0:
            continue
        to_Concatenate += ' '+line.strip()

outfile = args.jobs_dir+'ConcatenateGvcfs_chr'+args.chr+'.sh'
with open(outfile,'w') as out:
    new_template = template.replace('REPLACECHROMOSOME',args.chr)
    new_template = new_template.replace('REPLACEOUTDIR',args.outdir)
    new_template = new_template.replace('REPLACEINPUT', to_Concatenate)
    new_template = new_template.replace('REPLACEPROJECT',args.project)
    new_template = new_template.replace('REPLACEFINISHED','ConcatenateGvcfs_chr'+args.chr+'.sh.finished')
    out.write(new_template)
    print(outfile)

