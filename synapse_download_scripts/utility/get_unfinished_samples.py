import synapseclient
import getpass
import synapseutils
import glob

user = input("Synapse username:")
password = getpass.getpass('Synapse password:')
syn = synapseclient.Synapse()
syn.login(user,password)

samples_processed = {'featureCounts':set([]), 'kallisto':set([]), 'rMats':set([])}
studies = ['MayoTCX','MayoCBE','MSBB','ROSMAP']
for study in studies:
    for featureCount_file in glob.iglob('featureCounts/results/'+study+'/exon.countAll/*.txt.gz'):
        sample_name = featureCount_file.split('/')[-1].split('.exon.countAll')[0]
        samples_processed['featureCounts'].add(sample_name)

    for kallisto_file in glob.iglob('kallisto/results/'+study+'/*/abundance.tsv'):
        x = 0
        has_output = False
        with open(kallisto_file) as input_file:
            for line in input_file:
                x += 1
                if x >= 2:
                    has_output = True
                    break
        if not has_output:
            continue
        sample_name = kallisto_file.split('/')[-2]
        samples_processed['kallisto'].add(sample_name)
    
    for rMats_file in glob.iglob('rMats/rMats_results/'+study+'/*/A3SS.MATS.JCEC.txt'):
        x = 0
        has_output = False
        with open(rMats_file) as input_file:
            for line in input_file:
                x += 1
                if x >= 2:
                    has_output = True
                    break
        sample_name = rMats_file.split('/')[-2]
        if not has_output:
            print('rm '+sample_name+'.err')
            continue
        samples_processed['rMats'].add(sample_name)

def fix_name(name):
    name = name.split('.accepted')[0]
    name = name.split('Aligned')[0]
    name = name.split('.')[0]
    name = name.split('_resequenced')[0]
    return(name)

for program in samples_processed:
    with open('samples_not_processed.'+program+'.txt','w') as out:
        out.write('sample\tstudy\tsynapse_id\n')
        for file in syn.getChildren('syn8540822'):
            fixed_name = fix_name(file['name'])
            if fixed_name not in samples_processed[program]:
                out.write(fixed_name+'\tMSBB\t'+file['id']+'\n')
#                if program == 'featureCounts':
#                    synapseutils.syncFromSynapse(syn, file['id'], path = '/groups/umcg-biogen/tmp04/input/rawdata/AMP_AD/BAMs/leftover/')

        for file in syn.getChildren('syn8540821'):
            fixed_name = fix_name(file['name'])
            if fixed_name not in samples_processed[program]:
                out.write(fixed_name+'\tMayoCBE\t'+file['id']+'\n')
#                if program == 'featureCounts':
#                    synapseutils.syncFromSynapse(syn, file['id'], path = '/groups/umcg-biogen/tmp04/input/rawdata/AMP_AD/BAMs/leftover/')

        for file in syn.getChildren('syn8540820'):
            fixed_name = fix_name(file['name'])
            if fixed_name not in samples_processed[program]:
                out.write(fixed_name+'\tMayoTCX\t'+file['id']+'\n')
#                if program == 'featureCounts':
#                    synapseutils.syncFromSynapse(syn, file['id'], path = '/groups/umcg-biogen/tmp04/input/rawdata/AMP_AD/BAMs/leftover/')

        for file in syn.getChildren('syn8540863'):
            fixed_name = fix_name(file['name'])
            if fixed_name not in samples_processed[program]:
                out.write(fixed_name+'\tROSMAP\t'+file['id']+'\n')
#                if program == 'featureCounts':
#                    synapseutils.syncFromSynapse(syn, file['id'], path = '/groups/umcg-biogen/tmp04/input/rawdata/AMP_AD/BAMs/leftover/')

