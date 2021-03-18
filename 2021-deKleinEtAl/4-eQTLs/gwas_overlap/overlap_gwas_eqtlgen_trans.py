import glob
import gzip 

disease_per_snp = {}
qtl_per_disease = {}
diseases = []
for f in glob.glob('GWAS*'):
    disease = f.split('-')[-1].split('.')[0]
    qtl_per_disease[disease] = []
    diseases.append(disease)
    with open(f) as input_file:
        header = input_file.readline()
        for line in input_file:
            line = line.strip().split('\t')
            snp = line[2]
            if snp not in disease_per_snp:
                disease_per_snp[snp] = set([])
            disease_per_snp[snp].add(disease)

with gzip.open('/groups/umcg-wijmenga/tmp03/projects/eQTLGen/analysis/trans-eqtl/pcCorrected/output/v1013/2018-01-24-CombinedMetaAnalysis/standardization_for_collaborators/2018-09-04-trans-eQTLsFDR-CohortInfoRemoved.txt.gz') as input_file:
    header = input_file.readline().decode('utf-8')
    for disease in diseases:
        qtl_per_disease[disease].append(header)
    for line in input_file:
        line = line.decode('utf-8')
        snp = line.split('\t')[1]
        fdr = float(line.split('\t')[-1])
        if snp in disease_per_snp and fdr < 0.05:
            for disease in disease_per_snp[snp]:
                qtl_per_disease[disease].append(line)

for disease in qtl_per_disease:
    with open('eqtlgen_disease_overlap/eqtlGen_eQTLsFDR0.05_trans-'+disease+'.txt','w') as out:
        out.write(''.join(qtl_per_disease[disease]))

