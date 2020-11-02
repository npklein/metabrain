from termcolor import colored
from scipy.stats import t
from scipy.stats import norm
import math
import re
import sqlite3
import argparse
from beautifultable import BeautifulTable
from datetime import datetime

##### !! Don't change values below !! ####
# Number of genes included in the coreg table
# Used for calculating bonferonni corrected p-value
n_genes_all = 59056
n_genes_amygdala = 53147
n_genes_basal = 58805
n_genes_cerebellum = 58805
n_genes_cortex = 58999
n_genes_hippo = 58052
n_genes_hypo = 56667
n_genes_spinal = 57386

# n_tests = (number of genes * number of genes) / 2 because it is symetrical matrix
# so only half of tests of full matrix
n_tests = {'all':(n_genes_all*n_genes_all)/2,'amyg':(n_genes_amygdala*n_genes_amygdala)/2,
                'basal':(n_genes_basal*n_genes_basal)/2,'cerebellum':(n_genes_cerebellum*n_genes_cerebellum)/2,
                'cortex':(n_genes_cortex*n_genes_cortex)/2,'hippocampus':(n_genes_hippo*n_genes_hippo)/2,
                'hypothalamus':(n_genes_hypo*n_genes_hypo)/2,'spinal':(n_genes_spinal*n_genes_spinal)/2}

# Number of eigenvectors used to make the coreg table
# These are used to calculate the z-scores from the correlation
n_PCs = {'all':1000,'amyg':74, 'basal':200, 'cerebellum':300,
             'cortex':500, 'hippocampus':100, 'hypothalamus':200,
             'spinal':200}
##### !! Don't change values above !! ####




class GetCoreg():
    def __enter__(self):
        # Connection the the database
        self.conn = sqlite3.connect('MetaBrain.correlation.all-regions.db')
        self.c = self.conn.cursor()
        return self

    def __init__(self, db_location=None):
        self.db_location = db_location

    def get_coreg(self,gene1, gene2, regions=['all','amyg','basal','cerebellum',
                                         'cortex','hippocampus','hypothalamus',
                                         'spinal']):
        if not self.db_location:
            raise RutimeError('db location was not set')
        allowed_regions = set(['all','amyg','basal','cerebellum',
                                         'cortex','hippocampus','hypothalamus',
                                         'spinal'])
        for e in regions:
            if e not in allowed_regions:
                raise RuntimeError('regions included "'+e+'". Regions can only include: '+','.join(allowed_regions))
        gene_pair = (shortenENSG(gene1)+'_'+shortenENSG(gene2),)
        co_regions_names = 'cor_'+',cor_'.join(regions)
        if self.args.debug:
            before = datetime.now()
        self.c.execute("SELECT "+co_regions_names+" FROM correlations where gene_pair=?", 
                    gene_pair)
        if self.args.debug:
            diff = datetime.now()-before
            print('took',diff.total_seconds(),'seconds')

        row = self.c.fetchone()
        if not row:
            gene_pair_reverse = (shortenENSG(gene2)+'_'+shortenENSG(gene1),)
            self.c.execute("SELECT "+co_regions_names+" FROM correlations where gene_pair=?",
                      gene_pair_reverse)
            row = self.c.fetchone()
            if not row:
                raise RuntimeError('Gene pair '+gene1+' + '+gene2+' not in database')
        return(row)

    def __exit__(self, type, value, traceback):
        self.conn.close()

    def command_line(self):
        '''Use command line arguments to print coreg table
           This function runs if program is run from the command line
        '''
        parser = argparse.ArgumentParser(description='Select correlation between two genes.')
        parser.add_argument('gene1', help='1st gene in pair')
        parser.add_argument('gene2', help='2nd  gene in pair')
        parser.add_argument('db_location', help='database location')
        parser.add_argument('--regions', help='comma separated regions to select',
                            default='all,amyg,basal,cerebellum,'+
                                     'cortex,hippocampus,hypothalamus,'+
                                     'spinal')
        parser.add_argument('--debug', action='store_true')
        self.db_location = args.db_location
        args = parser.parse_args()
        self.args = args
        if type(args.regions) == str:
            args.regions = args.regions.split(',')
        coreg_scores, pvalues, pvalues_bonf, zscores = self.get_coreg_and_zScores(args.gene1, args.gene2, args.regions)

        print('Coreg between '+args.gene1 +' and '+args.gene2)

        table = BeautifulTable(precision=332)
        table.set_style(BeautifulTable.STYLE_COMPACT)
        for index, region in enumerate(args.regions):
            pval =  "{:.2e}".format(pvalues[index])
            pval_bonf  = "{:.2e}".format(pvalues_bonf[index])
            coreg = round(coreg_scores[index], 3)
            zscore = str(round(zscores[index], 3))[0:6]
            if pvalues[index] < 0.05:
                pval = colored(pval, "green")
            if pvalues_bonf[index] < 0.05:
               pval_bonf =  colored(pval_bonf,"green",attrs=["underline"])

            table.rows.append([region, coreg, zscore,
                              pval, pval_bonf])

        table.columns.header = [colored('region',attrs=['bold']),colored('coreg',attrs=['bold']),colored('z-score',attrs=['bold']),
                                colored('p-value',attrs=['bold']),colored('p-value bonf corrected',attrs=['bold'])]
        print(table)

    def get_coreg_and_zScores(self,gene1, gene2, regions=['all','amyg','basal','cerebellum',
                                         'cortex','hippocampus','hypothalamus',
                                         'spinal']):
        coreg_scores = self.get_coreg(gene1, gene2, regions)
        pvalues = []
        pvalues_bonf = []
        zscores = []
        for index, region in enumerate(regions):
            p_val, zScore = correlationToZ(coreg_scores[index], n_PCs[region])
            pvalues.append(p_val)
            zscores.append(zScore)
            bonf_cor = p_val*n_tests[region]
            if bonf_cor > 1:
                bonf_cor = 1
            pvalues_bonf.append(bonf_cor)

        return(coreg_scores, pvalues, pvalues_bonf, zscores)


def correlationToZ(correlation, nrSamples):
    '''Convert correlation to zscores using sample size'''
    try:
        t_val = correlation / math.sqrt((1.0 - correlation * correlation) / (nrSamples - 2))
    except (ValueError, ZeroDivisionError):
        print(correlation, nrSamples)
        raise
    p_value = 0.0
    z_score = 0.0
    if t_val < 0.0:
        # *2 because .cdf is two-tailed, but this is 1-tailed test
        p_value = t.cdf(t_val, nrSamples-2) * 2
        if p_value < 2.0e-323:
            p_value = 2.0e-323
        z_score = norm.ppf(p_value);
    else:
        # *2 because .cdf is two-tailed, but this is 1-tailed test
        p_value = t.cdf(-t_val, nrSamples-2) *2
        if p_value < 2.0e-323:
            p_value = 2.0e-323
    
        z_score = -norm.ppf(p_value)
    return p_value,z_score

def shortenENSG(string):
    '''Function to replace all leading
       zeros from a a given string
       to a number of replaced zeros
        e.g. 00023404 becomes 323404
    '''
    if string.startswith('ENSG'):
        string = string.replace('ENSG','')
    nzeros = len(string) - len(string.lstrip('0'))
    # Regex to remove leading
    # zeros from a string
    regex = "^0+(?!$)"
    # Replaces the matched
    # value with given string
    string = str(nzeros) + re.sub(regex, "", string)
    return(string)



if __name__ == "__main__":
    with  GetCoreg() as get_coreg:
        get_coreg.command_line()
