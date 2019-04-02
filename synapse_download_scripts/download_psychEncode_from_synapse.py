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
#files = synapseutils.syncFromSynapse(syn, 'syn4588490', path = 'RNAseq/EpiGABA/')
#files = synapseutils.syncFromSynapse(syn, 'syn4588489', path = 'metadata/EpiGABA/')
# UCLA-ASD
#files = synapseutils.syncFromSynapse(syn, 'syn4587614', path = 'metadata/UCLA_ASD/')
#files = synapseutils.syncFromSynapse(syn, 'syn4587615', path = 'RNAseq/UCLA_ASD/')
# CMC HBCC
#files = synapseutils.syncFromSynapse(syn, 'syn10254168', path = 'RNAseq/CMC_HBCC/bam/')


# Download metadata for Brainseq(LIBD_szControl), CMC, and CMC_HBCC to check if we already have all those samples
#results = syn.tableQuery('select * from syn8466658 where "study" = \'CMC\' AND "assay" = \'rnaSeq\'')
#CMC_meta_dir = 'metadata/CMC/'
#if not os.path.exists(CMC_meta_dir):
#    os.makedirs(CMC_meta_dir)
#with open(CMC_meta_dir+'/CMC_RNAseq_metadata.txt','w') as out:
#    for row in results:
#        row = [str(x) for x in row]
#        out.write('\t'.join(row)+'\n')

#results = syn.tableQuery('select * from syn8466658 where "study" = \'CMC_HBCC\' AND "assay" = \'rnaSeq\'')
#CMC_HBCC_meta_dir = 'metadata/CMC_HBCC'
#if not os.path.exists(CMC_HBCC_meta_dir):
#    os.makedirs(CMC_HBCC_meta_dir)
#with open(CMC_HBCC_meta_dir+'/CMC_HBCC_RNAseq_metadata.txt','w') as out:
#    for row in results:
#        row = [str(x) for x in row]
#        out.write('\t'.join(row)+'\n')

#results = syn.tableQuery('select * from syn8466658 where "study" = \'LIBD_szControl\' AND "assay" = \'rnaSeq\'')
#libd_meta_dir = 'metadata/LIBD_szControl/'
#if not os.path.exists(libd_meta_dir):
#    os.makedirs(libd_meta_dir)
#with open(libd_meta_dir+'/LIBD_szControl_RNAseq_metadata.txt','w') as out:
#    for row in results:
#        row = [str(x) for x in row]
#        out.write('\t'.join(row)+'\n')

# synapse get -q "SELECT * FROM syn8466658 WHERE ( ( \"study\" = 'CMC' OR \"study\" = 'LIBD_szControl' OR \"study\" = 'CMC_HBCC' ) AND ( \"assay\" = 'rnaSeq' ) )"
#files = synapseutils.syncFromSynapse(syn, 'syn4587614', path = 'metadata/LIBD_szControl/')

