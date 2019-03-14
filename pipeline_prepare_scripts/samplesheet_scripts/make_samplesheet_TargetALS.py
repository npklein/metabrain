# Parse TargetALS samplesheet to make molgenis-compute samplesheet
import os.path
import glob

parser = argparse.ArgumentParser(description='Make Molgenis Compute samplesheet for TargetALS.')
parser.add_argument('samplesheet', help='TargetALS samplesheet')
parser.add_argument('fastq_base_dir', help='path of base dir where subdirs with fastq files are in (will be recursively searched')
parser.add_argument('output_prefix',help='Prefix of output file')

args = parser.parse_args()

info_per_quote = {}
with open(args.samplesheet) as input_file:
    header = input_file.readline()
    for line in input_file:
        line = line.strip().split('\t')
        sample_name = line[2]
        quote = line[0]
        external_sample_id = line[1]
        external_subject_id = line[4]

        # have to add a lot of wildcards because the naming is not uniform and there are
        # variable parts in the name that can not be known before hand
        fastq1 = args.fastq_base_dir+'/'+quote+'/Project_*/Sample_'+external_sample_id.replace('_HRA_','-HRA-')
        fastq1 += '*/fastq/'+external_sample_id.replace('_HRA_','-HRA-')+'*'
        fastq2 = fastq1 + '.R2.fastq.gz'
        fastq1 += '.R1.fastq.gz'

        # search using wildcards because sample info does not include everything to build the full path
        # (for example the date can not be known from just the samplesheet)
        # Since fastq files get deleted, the number of the protocols in compute can change if samples not added
        # when the fastq file does not exist. Therefore, instead set fastq to NOT_FOUND if the fastq file can not be found
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

for quote in info_per_quote
    # don't add / after output_prefix because is part of name, not dir:
    with open(output_prefix+quote+'.txt', 'w') as out:
        out.write('internalId,project,sampleName,reads1FqGz,reads2FqGz,cram\n')
        for info in info_per_quote[quote]:
            out.write(info)
