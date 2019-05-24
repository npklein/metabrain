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
download_MAYOTCX = True
download_MAYOCBE = False
download_ROSMAP = False#True
download_all = False

# STAR MSBB
if download_MSBB:
    if not download_all:
        for file in syn.getChildren('syn8540822'):
#            if file['name'].startswith('hB'):
#                synapseutils.syncFromSynapse(syn, file['id'], path = 'BAMs/MSBB/')
            print(file['id'], file['name'])
        exit()
    else:
        # aligned BAM files MSBB
        files = synapseutils.syncFromSynapse(syn, 'syn8540822', path = 'BAMs/MSBB/')
        # STAR MSBB
        files = synapseutils.syncFromSynapse(syn, 'syn12104381', path = 'MSBB/STAR')
        # metadata MSBB
        files = synapseutils.syncFromSynapse(syn, 'syn7392158', path = 'metadata/')

if download_MAYOCBE:
    # STAR MAYO
    files = synapseutils.syncFromSynapse(syn, 'syn12104376', path = 'MAYO/STAR/')
    # aligned BAM files Mayo
    files = synapseutils.syncFromSynapse(syn, 'syn8540821', path = 'BAMs/MayoCBE/')

if download_MAYOTCX:
    files = synapseutils.syncFromSynapse(syn, 'syn8540820', path = 'BAMs/MayoTCX/')

    # metadata MAYO
    files = synapseutils.syncFromSynapse(syn, 'syn11384571', path = 'metadata/')
    files = synapseutils.syncFromSynapse(syn, 'syn5223705', path = 'metadata/')
    files = synapseutils.syncFromSynapse(syn, 'syn3817650', path = 'metadata/')


to_download = ['40_120416', '402_120503', '403_120503', '410_120503', '414_120503', '415_120503', '429_120507', '406_120503', '408_120503', '411_120503', '416_120503', '418_120507', '425_120507', '428_120507', '944_131107', '407_120503', '412_120503', '413_120503', '419_120507', '424_120507', '427_120507', '405_120503', '420_120507', '423_120507', '426_120507']
if download_ROSMAP:
    download_all = False
    # aligned bam files ROSMAP
    # only download certain files
    if not download_all:
        for file in syn.getChildren('syn8540863'):
            if file['name'].split('Aligned')[0] in to_download:
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

