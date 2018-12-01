# Parse TargetALS samplesheet to make molgenis-compute samplesheet
import sys
import os.path
import glob

if len(sys.argv) != 3:
    raise RuntimeError('Not correct number of arguments')

info_per_quote = {}
with open(sys.argv[1]) as input_file:
    header = input_file.readline()
    for line in input_file:
        line = line.strip().split('\t')
        sample_name = line[2]
        quote = line[0]
        external_sample_id = line[1]
        external_subject_id = line[4]

        # have to add a lot of wildcards because the naming is not uniform and there are
        # variable parts in the name that can not be known before hand
        baseDir = '/groups/umcg-biogen/tmp04/biogen/input/TargetALS/'
        fastq1 = baseDir+quote+'/Project_*/Sample_'+external_sample_id.replace('_HRA_','-HRA-')
        fastq1 += '*/fastq/'+external_sample_id.replace('_HRA_','-HRA-')+'*'
        fastq2 = fastq1 + '.R2.fastq.gz'
        fastq1 += '.R1.fastq.gz'

        # search using wildcards because sample info does not include everything to build the full path
        # (for example the date can not know)
        try:
            fastq1 = glob.glob(fastq1)[0]
        except IndexError:
            print(fastq1 + ' does not exist')
            fastq1 = "NOT_FOUND"
        try:
            fastq2 = glob.glob(fastq2)[0]
        except:
            print(fastq2+' does not exist')
            fastq2 = "NOT_FOUND"
#            raise RuntimeError(fastq1 + ' or '+fastq2+' does not exist')

        if quote not in info_per_quote:
            info_per_quote[quote] = []
        external_subject_id = external_subject_id.replace(' ','_')
        info_per_quote[quote].append(external_sample_id+',TargetALS,'+external_subject_id+','+fastq1+','+fastq2+','+'${cramFileDir}/${uniqueID}.cram\n')

for quote in info_per_quote:
    with open(sys.argv[2]+quote+'.txt', 'w') as out:
        out.write('internalId,project,sampleName,reads1FqGz,reads2FqGz,cram\n')
        for info in info_per_quote[quote]:
            out.write(info)
