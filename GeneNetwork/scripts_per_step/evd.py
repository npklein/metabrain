import sys
import os
import pandas as pd
import numpy as np
import argparse
from sklearn.decomposition import PCA
from pathlib import Path
from datetime import datetime

parser = argparse.ArgumentParser(description='Do PCA over correlation matrix')
parser.add_argument('corfile',help='Path to correlation file')
parser.add_argument('expressionfile',help='Path to expression file that was used to make the correlation file')
parser.add_argument('outdir',help='Path to the output directory')
parser.add_argument('--svd_solver',
                    help='svd solver to use (see https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html)',
                    default='auto')
args = parser.parse_args()

df = pd.read_csv(
    filepath_or_buffer=args.corfile,
    sep='\t',index_col=0)

Path(args.outdir).mkdir(parents=True, exist_ok=True)

df = pd.read_csv(
    filepath_or_buffer=args.expressionfile,
    sep='\t',index_col=0)


pca = PCA()
#pca = PCA(args.svd_solver)
now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")

print('Start PCA - '+dt_string)
sys.stdout.flush()
projected_data = pca.fit_transform(df.T)
eigenvalues = pca.explained_variance_
print('done')
sys.stdout.flush()
components = pca.components_
eigenvectors = pd.DataFrame(components.T)
pc_scores = pd.DataFrame(projected_data)
eigenvectors.columns = eigenvectors.columns + 1
eigenvectors = eigenvectors.add_prefix("PC")
eigenvectors.index = df.index
eigenvectors.index.name = datetime.now().strftime('%d/%m/%Y')
pc_scores.columns = pc_scores.columns + 1
pc_scores = pc_scores.add_prefix("PC")
pc_scores.index = df.columns
pc_scores.index.name = datetime.now().strftime('%d/%m/%Y')
eigenvectors.to_csv(os.path.join(args.outdir, "eigenvectors.txt"),sep='\t')
pc_scores.to_csv(os.path.join(args.outdir, "pc-scores.txt"),sep='\t')


#print('Calculate cronbach alpha')
#def CronbachAlpha(pc_scores, eigenvectors):
#    n_samples = eigenvectors.shape[0]
#    # Calculates Cronbach's alpha values for each component
    # Based on the cronbach alpha implementation in pca.cc of pca++
    # Only works if evd was calculated on correlation matrix
#    print("NOTE! Cronbach alpha calculation is assuming that evd is done on correlation matrix")
#    eigenvector_sums = eigenvectors.sum(axis=0)
#    pc_var = pc_scores.var(axis=0)
#    alphas = (n_samples / (n_samples - 1.0)) * (1.0 - eigenvector_sums.values / pc_var.values)
#    return(alphas)

eigenvalues = pd.DataFrame(eigenvalues, columns=['eigenvalues']).to_csv('eigenvalues.txt', index=False)

#cronbach = CronbachAlpha(pc_scores, eigenvectors)
#cronbach = pd.DataFrame(cronbach, columns=['cronbach']).to_csv('cronbach.txt', index=False)

