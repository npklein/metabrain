import os
import synapseclient
import synapseutils
import argparse

parser = argparse.ArgumentParser(description='Download RNAseq and genotypes of CMC.')
parser.add_argument('RNAseq_directory', help='Directory to download RNAseq data to')

args = parser.parse_args()

user = raw_input("Synapse username:")
password = getpass.getpass('Synapse password:')

syn = synapseclient.Synapse()
syn.login(user,password)

print('Sync Brainseq')
# RNAseq
files = synapseutils.syncFromSynapse(syn, 'syn8227833', path = 'RNAseq/')

