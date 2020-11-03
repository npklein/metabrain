
import gzip
terms = []
count_term = {}
term_type = {}
print('count go')
with open('term_order.txt','w') as out:
    out.write('term\torder\n')
    x = 1
    with gzip.open('2020-03-29-goa_human.gaf.gz','rt') as input_file:
        for index, line in enumerate(input_file):
            if line.startswith('!'):
                continue
            line = line.strip().split('\t')
            type = None
            if line[8] == 'F':
                type = 'GO:MF'
            elif line[8] == 'P':
                type = 'GO:BP'
            elif line[8] == 'C':
                type = 'GO:CC'
            else:
                raise RuntimeError(line[8])
            if line[4] not in count_term:
                count_term[line[4]] = 0
                terms.append(line[4])
                term_type[line[4]] = type
                out.write(line[4]+'\t'+str(x)+'\n')
                x += 1
            count_term[line[4]] += 1

    print('count kegg')
    with gzip.open('2020-03-28-c2.cp.kegg.v7.0.entrez.gmt.gz','rt') as input_file:
        for index, line in enumerate(input_file):
            line = line.split('\t')
            count_term[line[0]] = len(line)-2
            terms.append(line[0])
            term_type[line[0]] = 'KEGG'
            out.write(line[0]+'\t'+str(index)+'\n')
        

    print('count reactome')
    with gzip.open('2020-03-28-Ensembl2Reactome_All_Levels.txt.gz','rt') as input_file:
        x = 1
        for line in input_file:
            line = line.strip().split('\t')
            if line[1] not in count_term:
                count_term[line[1]] = 0
                terms.append(line[1])
                term_type[line[1]] = 'REAC'
                out.write(line[1]+'\t'+str(x)+'\n')
                x += 1
            count_term[line[1]] += 1


    print('count hpo')
    with gzip.open('2020-03-28-HPO-phenotype-to-genes.txt.gz','rt') as input_file:
        x = 1
        for line in input_file:
            if line.startswith('#'):
                continue
            line = line.strip().split('\t')
            if line[0] not in count_term:
                count_term[line[0]] = 0
                terms.append(line[0])
                out.write(line[0]+'\t'+str(x)+'\n')
                term_type[line[0]] = 'HP'
                x += 1
            count_term[line[0]] += 1


    with open('ngenes_per_term.txt','w') as out:
        out.write('term_id\tterm_size\tterm_type\n')
        for term in terms:
            out.write(term+'\t'+str(count_term[term])+'\t'+term_type[term]+'\n')




