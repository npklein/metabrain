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

args = parser.parse_args()

template = """#!/bin/bash
#SBATCH --job-name=GenotypeGvcf_chrREPLACECHROMOSOME
#SBATCH --output=GenotypeGvcf_chrREPLACECHROMOSOME.out
#SBATCH --error=GenotypeGvcf_chrREPLACECHROMOSOME.err
#SBATCH --time=23:59:59
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

$EBROOTGATK/gatk GenotypeGVCFs \\
    --reference REPLACEREFGENOME \\
    --dbsnp REPLACEDBSNP \\
    --output REPLACEOUTDIR/REPLACEOUTPUT \\
    --variant gendb://REPLACEINPREFIX.chrREPLACECHROMOSOME \\
    --stand-call-conf 10.0 \\
    -L chrREPLACECHROMOSOME \\
    -G StandardAnnotation

if [ "$?" -eq 0 ];
then
 echo "returncode: $?"; 
 
 cd 
 md5sum REPLACEOUTPUT > REPLACEOUTPUT.md5
 cd -
 echo "succes moving files";
else
 echo "returncode: $?";
 echo "fail";
fi

echo "## "$(date)" ##  $0 Done "I

"""

chromosomes = ['1','2','3','4','5','6','7',
        '8','9','10','11','12','13','14',
         '15','16','17','18','19','20',
         '21','22','X','Y','M']

for chr in chromosomes:
    outfile = args.jobs_dir+'GenotypeGvcf_chr'+chr+'.sh'
    with open(outfile,'w') as out:
        new_template = template.replace('REPLACECHROMOSOME',chr)
        new_template = new_template.replace('REPLACEOUTDIR',args.outdir)
        new_template = new_template.replace('REPLACEOUTPUT',args.project+'_chr'+chr+'.gg.vcf.gz')
        new_template = new_template.replace('REPLACEINPREFIX', args.genomicDB_dir+args.project)
        new_template = new_template.replace('REPLACEREFGENOME',args.ref_genome)
        new_template=new_template.replace('REPLACEDBSNP', args.dbsnp)
        out.write(new_template)
    print(outfile)

