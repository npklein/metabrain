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
download_MSBB = False
download_MAYO = False
download_ROSMAP = True

# STAR MSBB
if download_MSBB:
    files = synapseutils.syncFromSynapse(syn, 'syn12104381', path = 'MSBB/STAR')
    # metadata MSBB
    files = synapseutils.syncFromSynapse(syn, 'syn7392158', path = 'metadata/')
    # aligned BAM files MSBB
    files = synapseutils.syncFromSynapse(syn, 'syn8540822', path = 'BAMs/MSBB/')

if download_MAYO:
    # STAR MAYO
    files = synapseutils.syncFromSynapse(syn, 'syn12104376', path = 'MAYO/STAR/')
    # aligned BAM files Mayo
    files = synapseutils.syncFromSynapse(syn, 'syn8540821', path = 'BAMs/MayoCBE/')
    files = synapseutils.syncFromSynapse(syn, 'syn8540820', path = 'BAMs/MayoTCX/')

    # metadata MAYO
    files = synapseutils.syncFromSynapse(syn, 'syn11384571', path = 'metadata/')
    files = synapseutils.syncFromSynapse(syn, 'syn5223705', path = 'metadata/')
    files = synapseutils.syncFromSynapse(syn, 'syn3817650', path = 'metadata/')


if download_ROSMAP:
    download_all = False
    # aligned bam files ROSMAP
    # only download certain files
    if not download_all:
        for file in syn.getChildren('syn8540863'):
            if file['name'].startswith('4'):
                synapseutils.syncFromSynapse(syn, file['id'], path = 'BAMs/ROSMAP/')
                print(file['id'], file['name'])
    else:
        # STAR ROSMAP
        files = synapseutils.syncFromSynapse(syn, 'syn12104384', path = 'ROSMAP/STAR/')
        files = synapseutils.syncFromSynapse(syn, 'syn8540863', path = 'BAMs/ROSMAP/')


    # ROSMAP metadata
    files = synapseutils.syncFromSynapse(syn, 'syn3157322', path = 'metadata')
    files = synapseutils.syncFromSynapse(syn, 'syn11958660', path = 'metadata')
    files = synapseutils.syncFromSynapse(syn, 'syn11384589', path = 'metadata')

