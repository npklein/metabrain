
import pandas as pd; 
import argparse

parser = argparse.ArgumentParser(description='Compare eQTLgen to MetaBrain')
parser.add_argument('metaBrain_cis_eqtl',
                    help='MetaBrain cis eQTL data')
parser.add_argument('trait',
                    help='GWAS trait')
parser.add_argument('outfile',
                    help='outfile name')

data=pd.read_csv(args.metaBrain_cis_eqtl,sep='\t')

target_data=data.loc[data['GWAS_TRAIT'].astype(str).str.contains(args.trait),:]; 

target_data.to_csv(path_or_buf=args.outfile,index=False, sep='\t')

