# made by Omar El Garwany
# modified by Niek de Klein
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Compare eQTLgen to MetaBrain')
parser.add_argument('metaBrain_data', 
                    help='MetaBrain GWAS eQTL data')
parser.add_argument('eQTLgen_cis_data', 
                    help='cis-eQTLs from eQTL gen')
parser.add_argument('outfile', 
                    help='outfile name')

args = parser.parse_args()
print(args.accumulate(args.integers))

metabrain_data=pd.read_csv(args.metaBrain_data, sep='\t')
print('Loaded meta-brain data')
compiled_data=pd.DataFrame()

chunk_num=0
for eqtlgen_chunk in pd.read_csv(args.eQTLgen_cis_data, sep='\t', chunksize=1000000):
	merged_data=pd.DataFrame()
	merged_data=metabrain_data.merge(eqtlgen_chunk, on=['SNPName'], 
                                                    how='left').sort_values(by=['SNPName', 'ProbeName_x','ProbeName_y','DISEASE', 'CISTRANS', 'PC', 'FX']).drop_duplicates(subset=['SNPName', 'ProbeName_x','ProbeName_y', 'DISEASE', 'CISTRANS', 'PC', 'FX'])
	compiled_data=compiled_data.append(merged_data)
	chunk_num = chunk_num + 1
	print("Chunk", chunk_num, "processed.")

#print('Loaded eqtlgen data')
with open(args.outfile, 'w') as out:
    out.write(compiled_data.to_csv(sep='\t', index=False))
#print('Joined and written data')

