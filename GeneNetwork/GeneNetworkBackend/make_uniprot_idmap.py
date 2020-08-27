from datetime import datetime
import gzip
today = datetime.now().strftime("%Y-%m-%d")

outfile = today+'-uniprot_ensembl_map.txt.gz'
print('Start reading and writing, can take sometime')
with gzip.open('idmapping_selected.tab.gz','rt') as input_file, gzip.open(outfile,'wt') as out:
    out.write('uniprot_id\tensembl_id\n')
    for line in input_file:
        line = line.rstrip('\n').split('\t')
        ensembl_id = line[18]
        uniprot_id = line[0]
        if line[1].split('_')[1] != 'HUMAN':
            continue
        if len(ensembl_id.strip()) > 0:
            out.write(uniprot_id+'\t'+ensembl_id+'\n')

print('output written to '+outfile)
