import os
import synapseclient
import synapseutils
print(synapseclient.__file__)

parser = argparse.ArgumentParser(description='Download RNAseq and genotypes of CMC.')
parser.add_argument('username', help='synapse username')
parser.add_argument('password', help='synapse password')
parser.add_argument('RNAseq_directory', help='Directory to download RNAseq data to')
parser.add_argument('Genotype_directory', help='Directory to download genotypes to')

args = parser.parse_args()

print('Sync AMP-AD')

# STAR MSBB
files = synapseutils.syncFromSynapse(syn, 'syn12104381', path = 'MSBB/STAR')#, exclude='.bam')
# STAR MAYO
files = synapseutils.syncFromSynapse(syn, 'syn12104376', path = 'MAYO/STAR/')#, exclude='.bam')
# STAR ROSMAP
files = synapseutils.syncFromSynapse(syn, 'syn12104384', path = 'ROSMAP/STAR/')#, exclude='.bam')

