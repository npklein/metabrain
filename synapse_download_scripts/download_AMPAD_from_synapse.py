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
# aligned BAM files MSBB
files = synapseutils.syncFromSynapse(syn, 'syn8540822', path = 'BAMs/MSBB/')#, exclude='.bam')
exit()
# STAR MAYO
#files = synapseutils.syncFromSynapse(syn, 'syn12104376', path = 'MAYO/STAR/')#, exclude='.bam')
# aligned BAM files Mayo
files = synapseutils.syncFromSynapse(syn, 'syn8540821', path = 'BAMs/MayoCBE/')#, exclude='.bam')
files = synapseutils.syncFromSynapse(syn, 'syn8540820', path = 'BAMs/MayoTCX/')#, exclude='.bam')

# metadata MAYO
files = synapseutils.syncFromSynapse(syn, 'syn11384571', path = 'metadata/')#, exclude='.bam')
files = synapseutils.syncFromSynapse(syn, 'syn5223705', path = 'metadata/')#, exclude='.bam')
files = synapseutils.syncFromSynapse(syn, 'syn3817650', path = 'metadata/')#, exclude='.bam')

# STAR ROSMAP
#files = synapseutils.syncFromSynapse(syn, 'syn12104384', path = 'ROSMAP/STAR/')#, exclude='.bam')
# aligned bam files ROSMAP
files = synapseutils.syncFromSynapse(syn, 'syn8540863', path = 'BAMs/ROSMAP/')#, exclude='.bam')
# ROSMAP metadata
files = synapseutils.syncFromSynapse(syn, 'syn3157322', path = 'metadata')#, exclude='.bam')
files = synapseutils.syncFromSynapse(syn, 'syn11958660', path = 'metadata')#, exclude='.bam')
files = synapseutils.syncFromSynapse(syn, 'syn11384589', path = 'metadata')#, exclude='.bam')

