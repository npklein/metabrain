import math
import glob  
import argparse

parser = argparse.ArgumentParser(description='Create slurm jobs to call genotype on merged gVCFs (all batches within the directory)')
parser.add_argument('merged_vcf_dir', help='Directory that contains the merged gVCFs')
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
#SBATCH --qos=leftover
#SBATCH --time=3-23:59:59
#SBATCH --cpus-per-task 10
#SBATCH --mem 100gb
#SBATCH --nodes 1
#SBATCH --export=NONE
#SBATCH --get-user-env=L


ENVIRONMENT_DIR="."
set -e
set -u

echo "## "$(date)" Start $0"

#Generate input files, according to number of batches
for i in {0..REPLACENUMBEROFBATCHES}
do
    inputs+=" --variant REPLACEINPREFIX.batch${i}_chrREPLACECHROMOSOME.g.vcf.gz"
done


#Load gatk module
module load GATK/3.8-0-Java-1.8.0_121
module list

mkdir -p REPLACEOUTDIR

if java -Xmx95g -XX:ParallelGCThreads=6 -Djava.io.tmpdir=${TMPDIR} -jar $EBROOTGATK/GenomeAnalysisTK.jar \\
 -T GenotypeGVCFs \\
 -R REPLACEREFGENOME \\
 --dbsnp REPLACEDBSNP \\
 -o REPLACEOUTDIR/REPLACEOUTPUT \\
 $inputs \\
 -stand_call_conf 10.0 \\
 -L chrREPLACECHROMOSOME
then
 echo "returncode: $?"; 
 
 cd REPLACEOUTDIR/
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

max_batch = 0
for f in glob.glob(args.merged_vcf_dir+'/*vcf.gz'):
    # NOTE: only works if the file naming is done in same way we do it
    batch = f.split('.')[1].split('_')[0]
    n = int(batch.replace('batch',''))
    if n > max_batch:
           max_batch = n

for chr in chromosomes:
    outfile = args.jobs_dir+'GenotypeGvcf_chr'+chr+'.sh'
    with open(outfile,'w') as out:
        new_template = template.replace('REPLACENUMBEROFBATCHES',str(max_batch))
        new_template = new_template.replace('REPLACECHROMOSOME',chr)
        new_template = new_template.replace('REPLACEOUTDIR',args.outdir)
        new_template = new_template.replace('REPLACEOUTPUT',args.project+'_chr'+chr+'.gg.vcf.gz')
        new_template = new_template.replace('REPLACEINPREFIX', args.merged_vcf_dir+'/'+args.project)
        new_template = new_template.replace('REPLACEREFGENOME',args.ref_genome)
        new_template=new_template.replace('REPLACEDBSNP', args.dbsnp)
        out.write(new_template)
    print(outfile)

