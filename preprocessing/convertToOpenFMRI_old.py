import os
import numpy as np
import pandas as pd
subList = ["LMVPA001", "LMVPA002", "LMVPA003", "LMVPA005",
           "LMVPA006", "LMVPA007", "LMVPA008", "LMVPA009",
           "LMVPA010", "LMVPA011", "LMVPA013", "LMVPA014",
           "LMVPA015", "LMVPA016", "LMVPA017", "LMVPA018", "LMVPA019", "LMVPA020"]
subData = {"LMVPA001": ["Run1", "Run2", "Run3", "Run4"],
           "LMVPA002": ["Run1", "Run2", "Run3", "Run4"],
           "LMVPA003": ["Run1", "Run2", "Run3", "Run4", "Run5", "Run6", "Run7", "Run8"],
           "LMVPA005": ["Run1", "Run2", "Run3", "Run4", "Run5", "Run6", "Run7", "Run8"],
           "LMVPA006": ["Run1", "Run2", "Run3", "Run4", "Run5", "Run6", "Run7", "Run8"],
           "LMVPA007": ["Run1", "Run2", "Run3", "Run4", "Run5", "Run6", "Run7", "Run8"],
           "LMVPA008": ["Run1", "Run2", "Run3", "Run4", "Run5", "Run6", "Run7", "Run8"],
           "LMVPA009": ["Run1", "Run2", "Run3", "Run4", "Run5", "Run6", "Run7", "Run8"],
           "LMVPA010": ["Run1", "Run2", "Run3", "Run4", "Run5", "Run6", "Run7", "Run8"],
           "LMVPA011": ["Run1", "Run2", "Run3", "Run4", "Run5", "Run6", "Run7", "Run8"],
           "LMVPA013": ["Run1", "Run2", "Run3", "Run4", "Run5", "Run6", "Run7", "Run8"],
           "LMVPA014": ["Run1", "Run2", "Run3", "Run4", "Run5", "Run6", "Run7", "Run8"],
           "LMVPA015": ["Run1", "Run2", "Run3", "Run4", "Run5", "Run6", "Run7", "Run8"],
           "LMVPA016": ["Run1", "Run2", "Run3", "Run4", "Run5", "Run6", "Run7", "Run8"],
           "LMVPA017": ["Run1", "Run2", "Run3", "Run4", "Run5", "Run6", "Run7", "Run8"],
           "LMVPA018": ["Run1", "Run2", "Run3", "Run4", "Run5", "Run6", "Run7", "Run8"],
           "LMVPA019": ["Run1", "Run2", "Run3", "Run4", "Run5", "Run6", "Run7", "Run8"]
           # "LMVPA020": ["Run1", "Run2", "Run3", "Run4", "Run5", "Run6", "Run7", "Run8"]
           }
# subData = {"LMVPA001": ["Run1", "Run2", "Run3", "Run4"]}

# projectDir = os.path.join(rootDir, 'CloudStation', 'scratchpad',

def writefile(sub, run):
    # screw it we start over
    rootDir = "D:\\"
    projectDir = os.path.join(rootDir, 'fmri', 'LanguageMVPA')
    t = (np.loadtxt(os.path.join(projectDir, 'regressor', sub, 'orig', sub + '_' + run + '.txt')))
    cnames = open(os.path.join(projectDir, 'regressor', sub, 'orig', sub + '_' + run + '_conds.txt'), "r")
    lines = cnames.readlines()
    ofile = open(os.path.join(projectDir, 'data', sub, 'func', sub + '_' + run + '.tsv'), "w")
    ofile.write('onset\tduration\ttrial_type\tstim\tverb\tanim\tAP\tCR\tsyntax\tprobe\n')
    for i, c in enumerate(lines):
        cid = c.split('_')[0].split(':')[1] + '_' + c.split('_')[1] + '_' + c.split('_')[3] \
               + '_' + c.split('_')[4] + '_' + c.split('_')[5].split(',')[0]
        ofile.write(str(t[i, 0]) + '\t' + str(t[i, 1]) + '\t' + cid + '\t' + addcondid(cid) + '\n')
    cnames.close()
    p = (np.loadtxt(os.path.join(projectDir, 'regressor', sub, 'orig', sub + '_' + run + '_probe.txt')))
    for i in np.arange(len(p)):
        ofile.write(str(p[i, 0]) + '\t' + str(p[i, 1]) + '\t' + 'probe' + '\t' + addcondid('probe') + '\n')
    ofile.close()

#take in string and convert to stuff
def addcondid(st):
    if st == 'probe':
        return '0\t0\t0\t0\t0\t0\t1'
    else:
        stype = np.array(['s', 'l'])
        verbs = np.array(['Touch', 'Light', 'Hit', 'Crush', 'Kiss', 'Stretch', 'Kick', 'Console'])
        syntax = np.array(['Act_Can', 'Act_Rel', 'Pass_Can', 'Pass_Rel'])
        ap = np.array(['Act', 'Pass'])
        cr = np.array(['Can', 'Rel'])
        anim = np.array(['Inanim', 'Anim'])

        return str(np.where(stype==st.split('_')[0])[0][0]+1) + '\t' + \
                str(np.where(verbs==st.split('_')[1])[0][0]+1) + '\t' + \
                str(np.where(anim==st.split('_')[2])[0][0]+1) + '\t' + \
                str(np.where(ap==st.split('_')[3])[0][0]+1) + '\t' + \
                str(np.where(cr==st.split('_')[4])[0][0]+1) + '\t' + \
                str(np.where(syntax==st.split('_')[3]+ '_' + st.split('_')[4])[0][0]+1) + '\t0'


# run the damn thing
for subject in subData.keys():
     [writefile(subject, r) for r in subData[subject]]