import os
import synapseclient
import synapseutils
import argparse

parser = argparse.ArgumentParser(description='Download RNAseq and genotypes of CMC.')
parser.add_argument('RNAseq_directory', help='Directory to download RNAseq data to')
parser.add_argument('Genotype_directory', help='Directory to download genotypes to')

args = parser.parse_args()


user = raw_input("Synapse username:")
password = getpass.getpass('Synapse password:')

syn = synapseclient.Synapse()
syn.login(user,password)

print('Sync CMC')

# RNAseq
if not os.path.exists(args.RNAseq_directory):
    os.makedirs(args.RNAseq_directory)
files = synapseutils.syncFromSynapse(syn, 'syn3280440', path = args.RNAseq_directory)
# Genotypes
if not os.path.exists(args.Genotype_directory):
    os.makedirs(args.Genotype_directory)
files = synapseutils.syncFromSynapse(syn, 'syn3275211', path = args.Genotype_directory)

