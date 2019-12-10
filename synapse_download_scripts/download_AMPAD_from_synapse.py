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
download_MSBB = True
download_MAYOTCX = True
download_MAYOCBE = True
download_ROSMAP = True
download_all = True

# STAR MSBB
if download_MSBB:
    if not download_all:
        for file in syn.getChildren('syn8540822'):
            if file['name'].startswith('hB'):
                synapseutils.syncFromSynapse(syn, file['id'], path = 'BAMs/MSBB/')
                print(file['id'], file['name'])
    else:
        # fastq files
        files = synapseutils.syncFromSynapse(syn, 'syn8612191', path = 'fastq/MSBB/')

        # aligned BAM files MSBB
#        files = synapseutils.syncFromSynapse(syn, 'syn8540822', path = 'BAMs/MSBB/')
        # STAR MSBB
#        files = synapseutils.syncFromSynapse(syn, 'syn12104381', path = 'MSBB/STAR')
        # metadata MSBB
        files = synapseutils.syncFromSynapse(syn, 'syn7392158', path = 'metadata/')

if download_MAYOCBE:
    # fastq
    files = synapseutils.syncFromSynapse(syn, 'syn8612213', path = 'fastq/MayoCBE/')


    # STAR MAYO
#    files = synapseutils.syncFromSynapse(syn, 'syn12104376', path = 'MAYO/STAR/')
    # aligned BAM files Mayo
#    files = synapseutils.syncFromSynapse(syn, 'syn8540821', path = 'BAMs/MayoCBE/')

if download_MAYOTCX:
#    files = synapseutils.syncFromSynapse(syn, 'syn8540820', path = 'BAMs/MayoTCX/')
    #fastq
    files = synapseutils.syncFromSynapse(syn, 'syn8612203', path = 'fastq/MayoTCX/')

    # metadata MAYO
    files = synapseutils.syncFromSynapse(syn, 'syn11384571', path = 'metadata/')
    files = synapseutils.syncFromSynapse(syn, 'syn5223705', path = 'metadata/')
    files = synapseutils.syncFromSynapse(syn, 'syn3817650', path = 'metadata/')


to_download = ['74_120417']
if download_ROSMAP:
    # aligned bam files ROSMAP
    # only download certain files
    if not download_all:
        for file in syn.getChildren('syn8540863'):
            if '120417' in file['name']:
                print(file['name'])
            if file['name'].split('Aligned')[0] in to_download:
                synapseutils.syncFromSynapse(syn, file['id'], path = 'BAMs/ROSMAP/')
                print(file['id'], file['name'])
    else:
        pass
        # STAR ROSMAP
#        files = synapseutils.syncFromSynapse(syn, 'syn12104384', path = 'ROSMAP/STAR/')
#        files = synapseutils.syncFromSynapse(syn, 'syn8540863', path = 'BAMs/ROSMAP/')

    #fastq
    files = synapseutils.syncFromSynapse(syn, 'syn8612097', path = 'fastq/ROSMAP/')

    # ROSMAP metadata
    files = synapseutils.syncFromSynapse(syn, 'syn3157322', path = 'metadata')
    files = synapseutils.syncFromSynapse(syn, 'syn11958660', path = 'metadata')
    files = synapseutils.syncFromSynapse(syn, 'syn11384589', path = 'metadata')

