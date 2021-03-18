# Download eQTLs from Sieberts et al. paper (https://www.biorxiv.org/content/biorxiv/early/2019/05/17/638544.full.pdf)
import synapseclient
import getpass
import synapseutils
from datetime import datetime

today = datetime.today().strftime('%Y-%m-%d')



user = input("Synapse username:")
password = getpass.getpass('Synapse password:')

syn = synapseclient.Synapse()
syn.login(user,password)


# ROSMAP eQTL
synapseutils.syncFromSynapse(syn, 'syn16984409', path = today+'-Sieberts-eQTLs/')

# Mayo TCX eQTL
synapseutils.syncFromSynapse(syn, 'syn16984410', path = today+'-Sieberts-eQTLs/')

# Mayo CER eQTL
synapseutils.syncFromSynapse(syn, 'syn16984411', path = today+'-Sieberts-eQTLs/')

# meta-analysis eQTL
synapseutils.syncFromSynapse(syn, 'syn16984815', path = today+'-Sieberts-eQTLs/')

