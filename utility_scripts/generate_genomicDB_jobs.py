
import math
import glob  
import argparse

parser = argparse.ArgumentParser(description='Create slurm jobs to make genomicDB from multiple gVCF files')
parser.add_argument('list_of_vcfs', help='.txt file with list of path of VCFs to merge (only for chr15, this will be replaced for other chr numbers)')
parser.add_argument('jobs_dir', help='directory to write jobs to')
parser.add_argument('outdir', help='directory to write merged gVCFs to')
parser.add_argument('project', help='Name of the project (will be prepended to the output file)')

args = parser.parse_args()

# Below is the job template that will be used
template = """#!/bin/bash
#SBATCH --job-name=genomicDB_chrREPLACECHROMOSOME
#SBATCH --output=genomicDB_chrREPLACECHROMOSOME.out
#SBATCH --error=genomicDB_chrREPLACECHROMOSOME.err
#SBATCH --time=2-23:59:59
#SBATCH --cpus-per-task 1
#SBATCH --mem 60gb
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

echo "NOTE: it crashes if the genomicDB already exists, so deleting it now first"
rm -rf REPLACEOUTDIR/REPLACEPROJECT.chrREPLACECHROMOSOME
$EBROOTGATK/gatk GenomicsDBImport \\
    --genomicsdb-workspace-path REPLACEOUTDIR/REPLACEPROJECT.chrREPLACECHROMOSOME \\
    --intervals /data/umcg-ndeklein/apps/data/storage.googleapis.com/gatk-test-data/intervals/hg38.even.handcurated.20k.chrREPLACECHROMOSOME.intervals \\
    REPLACEINPUT

if [ "$?" -eq 0 ];
then
  cd REPLACEOUTDIR/REPLACEPROJECT.chrREPLACECHROMOSOME
  for f in *;
  do
      [[ -e $f ]] || continue
      if [ -f "$f" ]; then
          md5sum $f > ${f}.md5
      fi
  done
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

input = {'1':'','2':'','3':'','4':'','5':'','6':'','7':'',
         '8':'','9':'','10':'','11':'','12':'','13':'','14':'',
         '15':'','16':'','17':'','18':'','19':'','20':'',
         '21':'','22':'','X':'','Y':'','M':''}


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
        input[str(chr)] += ' --variant '+line.strip().replace('chr15','chr'+str(chr))

for chr in input:
    outfile = args.jobs_dir+'genomicDB_chr'+chr+'.sh'
    with open(outfile,'w') as out:
        new_template = template.replace('REPLACECHROMOSOME',chr)
        new_template = new_template.replace('REPLACEOUTDIR',args.outdir)
        new_template = new_template.replace('REPLACEINPUT', input[chr])
        new_template = new_template.replace('REPLACEPROJECT',args.project)
        new_template = new_template.replace('REPLACEFINISHED','genomicDB_chr'+chr+'.sh.finished')
        out.write(new_template)
    print(outfile)

