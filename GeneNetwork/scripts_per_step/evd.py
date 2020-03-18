
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

print('Start PCA - '+dt_string, flush=True)
projected_data = pca.fit_transform(df.T)
eigenvalues = pca.explained_variance_
print('done', flush=True)
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


print('Calculate cronbach alpha')
def CronbachAlpha(itemscores):
    itemscores = np.asarray(itemscores)
    itemvars = itemscores.var(ddof=1)
    tscores = itemscores.sum()
    nitems = itemscores.size
    print(tscores)
    return nitems / (nitems-1.) * (1 - itemvars.sum() / tscores.var(ddof=1))

eigenvalues = pd.DataFrame(eigenvalues, columns=['eigenvalues']).to_csv('eigenvalues.csv')

for col in pc_scores.columns:
    cronbach = CronbachAlpha(pc_scores[pc_scores.columns[0]])
    print(cronbach)
cronbach = pd.DataFrame(cronbach, columns=['cronbach']).to_csv('cronbach.csv')

