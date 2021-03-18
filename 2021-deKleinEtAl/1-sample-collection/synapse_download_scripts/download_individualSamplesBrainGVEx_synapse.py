# Some samples from BrainGVEx had 0 counts after quantification, have to redo. This should be a one-time-thing so just using hard-coded paths

import os
import synapseclient
import synapseutils
import argparse
import getpass

user = input("Synapse username:")
password = getpass.getpass('Synapse password:')

syn = synapseclient.Synapse()
syn.login(user,password)

samples_to_redownload = set([])
with open('/groups/umcg-biogen/tmp04/input/rawdata/psychEncode/RNAseq/BrainGVEX/samples_with_0_counts.txt') as input_file:
    input_file.readline()
    for line in input_file:
        line = line.strip()
        samples_to_redownload.add(line)

syn_per_file = {}
with open('/groups/umcg-biogen/tmp04/input/rawdata/psychEncode/metadata/CMC_HBCC/CMC_HBCC_RNAseq_metadata.txt') as input_file:
    for line in input_file:
        line = line.strip().split('\t')
        if line[1].endswith('sortedByCoord.out.bam'):
            syn_per_file[line[1]] = line[0]

with open('/groups/umcg-biogen/tmp04/input/rawdata/psychEncode/metadata/BrainGVEX/BrainGVEx_RNA_table.txt') as input_file:
    for line in input_file:
        line = line.strip().replace('"','').split('\t')
        if line[3] in samples_to_redownload:
            files = synapseutils.syncFromSynapse(syn, line[0], path = 'redownload/')
