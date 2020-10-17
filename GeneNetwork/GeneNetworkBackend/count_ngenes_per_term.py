
import gzip
terms = []
count_term = {}
print('count go')
with gzip.open('2020-03-29-goa_human.gaf.gz','rt') as input_file:
    for line in input_file:
        if line.startswith('!'):
            continue
        line = line.strip().split('\t')
        if line[4] not in count_term:
            count_term[line[4]] = 0
            terms.append(line[4])
        count_term[line[4]] += 1

print('count kegg')
with gzip.open('2020-03-28-c2.cp.kegg.v7.0.entrez.gmt.gz','rt') as input_file:
    for line in input_file:
        line = line.split('\t')
        count_term[line[0]] = len(line)-2
        terms.append(line[0])

print('count reactome')
with gzip.open('2020-03-28-Ensembl2Reactome_All_Levels.txt.gz','rt') as input_file:
    for line in input_file:
        line = line.strip().split('\t')
        if line[1] not in count_term:
            count_term[line[1]] = 0
            terms.append(line[1])
        count_term[line[1]] += 1


print('count hpo')
with gzip.open('2020-03-28-HPO-phenotype-to-genes.txt.gz','rt') as input_file:
    for line in input_file:
        if line.startswith('#'):
            continue
        line = line.strip().split('\t')
        if line[0] not in count_term:
            count_term[line[0]] = 0
            terms.append(line[0])
        count_term[line[0]] += 1


with open('ngenes_per_term.txt','w') as out:
    out.write('term\tterm_size\n')
    for term in terms:
        out.write(term+'\t'+str(count_term[term])+'\n')
