import synapseclient
import argparse
import getpass
import synapseutils

parser = argparse.ArgumentParser(description='Download RNAseq and genotypes of CMC.')

args = parser.parse_args()


user = input("Synapse username:")
password = getpass.getpass('Synapse password:')

syn = synapseclient.Synapse()
syn.login(user,password)


print('Sync AMP-AD')

# STAR MSBB
# files = synapseutils.syncFromSynapse(syn, 'syn12104381', path = 'MSBB/STAR')#, exclude='.bam')
# metadata MSBB
files = synapseutils.syncFromSynapse(syn, 'syn7392158', path = 'metadata/')#, exclude='.bam')

# STAR MAYO
#files = synapseutils.syncFromSynapse(syn, 'syn12104376', path = 'MAYO/STAR/')#, exclude='.bam')
# metadata MAYO
files = synapseutils.syncFromSynapse(syn, 'syn11384571', path = 'metadata/')#, exclude='.bam')
files = synapseutils.syncFromSynapse(syn, 'syn5223705', path = 'metadata/')#, exclude='.bam')
files = synapseutils.syncFromSynapse(syn, 'syn3817650', path = 'metadata/')#, exclude='.bam')

# STAR ROSMAP
#files = synapseutils.syncFromSynapse(syn, 'syn12104384', path = 'ROSMAP/STAR/')#, exclude='.bam')
# ROSMAP metadata
files = synapseutils.syncFromSynapse(syn, 'syn3157322', path = 'metadata')#, exclude='.bam')
files = synapseutils.syncFromSynapse(syn, 'syn11958660', path = 'metadata')#, exclude='.bam')
files = synapseutils.syncFromSynapse(syn, 'syn11384589', path = 'metadata')#, exclude='.bam')

