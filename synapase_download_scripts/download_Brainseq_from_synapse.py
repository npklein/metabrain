import os
import synapseclient
import synapseutils
print(synapseclient.__file__)

parser = argparse.ArgumentParser(description='Download RNAseq and genotypes of CMC.')
parser.add_argument('RNAseq_directory', help='Directory to download RNAseq data to')
parser.add_argument('Genotype_directory', help='Directory to download genotypes to')

args = parser.parse_args()

user = raw_input("Synapse username:")
password = getpass.getpass('Synapse password:')

syn = synapseclient.Synapse()
syn.login(user,password)

print('Sync Brainseq')
# RNAseq
files = synapseutils.syncFromSynapse(syn, 'syn8227833', path = 'RNAseq/')

