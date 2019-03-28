import os
import synapseclient
import synapseutils
import argparse
import getpass

parser = argparse.ArgumentParser(description='Download RNAseq and genotypes of CMC.')
parser.add_argument('RNAseq_directory', help='Directory to download RNAseq data to')

args = parser.parse_args()

user = input("Synapse username:")
password = getpass.getpass('Synapse password:')

syn = synapseclient.Synapse()
syn.login(user,password)

print('Sync psychEncode')
# BipSeq
#files = synapseutils.syncFromSynapse(syn, 'syn8403872', path = 'RNAseq/BipSeq/')
#files = synapseutils.syncFromSynapse(syn, 'syn5845180', path = 'metadata/BipSeq/')
# BrainGVEX
#files = synapseutils.syncFromSynapse(syn, 'syn7062404', path = 'RNAseq/BrainGVEX/')
#files = synapseutils.syncFromSynapse(syn, 'syn3270014', path = 'metadata/BrainGVEX/')
# EpiGABA
files = synapseutils.syncFromSynapse(syn, 'syn4588490', path = 'RNAseq/EpiGABA/')
files = synapseutils.syncFromSynapse(syn, 'syn4588489', path = 'metadata/EpiGABA/')

