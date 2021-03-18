import glob

with open('copy_resources.sh','w') as out:
    for f in glob.glob('/groups/umcg-biogen/tmp03/input/rawdata/ucl-upload-biogen//Public_RNA-seq_*/parameter_files/parameters.csv'):
        with open(f) as input_file:
            header = input_file.readline()
            res_dir = ''
            for line in input_file:
                line = line.strip().split(',')
                if 'resDir' in ''.join(line):
                    if line[0] == 'resDir':
                        res_dir = line[1]
                    else:
                        line[1] = line[1].replace('${resDir}',res_dir)
                        # make the target directory
                        out.write('rsync -vP –rsync-path="mkdir -p /scratch/umcg-ndeklein/'+'/'.join(line[1].lstrip('/').split('/')[:-1])+' && rsync" '+line[1])
                        out.write(' umcg-ndeklein@peregrine.hpc.rug.nl:/scratch/umcg-ndeklein/'+'/'.join(line[1].lstrip('/').split('/')[:-1])+'/\n')
                elif '/apps/data' in ' '.join(line):
                    out.write('rsync -vPR –rsync-path="mkdir -p /scratch/umcg-ndeklein/'+'_'.join(line[1].lstrip('/').split('/')[:-1])+' && rsync" '+line[1])
                    out.write(line[1]+' umcg-ndeklein@peregrine.hpc.rug.nl:/scratch/umcg-ndeklein/'+'/'.join(line[1].lstrip('/').split('/')[:-1])+'/\n')
