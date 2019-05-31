import pandas as pd;
import argparse

parser = argparse.ArgumentParser(description='Compare eQTLgen to MetaBrain')
parser.add_argument('eQTL_file',
                    help='File with eQTL data')
parser.add_argument('gwas_file',
                    help='File with gwas data')
parser.add_argument('target_trait',
                    help='File with gwas data')
parser.add_argument('outfile',
                    help='outfile name')
args = parser.parse_args()

eqtl_data=pd.read_csv(args.eQTL_file, sep='\t');
gwas_data=pd.read_csv(args.gwas_file,header=None,sep='\t');
annotated_data=eqtl_data.merge(gwas_data, how='left', left_on=['SNPName'], right_on=[0]);
gwas_target_trait_num=gwas_data.loc[gwas_data[1].astype(str).str.contains(args.target_trait),:].shape[0];
gwas_total_num=gwas_data.shape[0];
annotated_eqtl_target_num=annotated_data.loc[annotated_data[1].astype(str).str.contains(args.target_trait),:].shape[0];
eqtl_total_num=eqtl_data.shape[0];
annotated_eqtl_num=annotated_data.loc[annotated_data[1].astype(str).str.contains(args.target_trait),:].shape[0];

with open(args.outfile,'w') as out:
    out.write(annotated_data.rename(index=str, columns={0: 'GWAS_SNP', 1: 'GWAS_TRAIT'}).to_csv(sep='\t', index=False,))

