# need to edit probe file to reflect that nothing actually happened in the probes
import os
import numpy as np
from collections import defaultdict


def writesubrun(f):
    sub = 'LMVPA' + str(f.split('_')[0])
    run = 'Run' + str(int(f.split('_')[1]))
    subData[sub].append(run)
    fullfile = np.loadtxt(os.path.join(projectDir, 'from_scanner', f), dtype=str, skiprows=1)
    ofile = open(os.path.join(projectDir, 'data', sub, 'func', sub + '_' + run + '.tsv'), "w")
    ofile.write('onset\tduration\ttrial_type\tstim\tverb\tanim\tAP\tCR\tsyntax\tprobe\n')
    for line in fullfile:
        resp = line.split(',')
        ofile.write(resp[1] + '\t' + resp[2] + '\t' + resp[0] + '\t' + addcondid(resp[0]) + '\n')
        if (resp[-1] == '') or (resp[-1] == 'None') or (resp[-1] == '0') or (resp[-1] == 0):
            ofile.write(resp[3] + '\t' + resp[4] + '\t' + 'probe' + '\t' + addcondid(None) + '\n')
        else:
            ofile.write(resp[3] + '\t' + resp[4] + '\t' + 'probe' + '\t' + addcondid('probe') + '\n')
    ofile.close()


def addcondid(st):
    if st is not None:
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
                    str(np.where(anim==st.split('_')[3])[0][0]+1) + '\t' + \
                    str(np.where(ap==st.split('_')[4])[0][0]+1) + '\t' + \
                    str(np.where(cr==st.split('_')[5])[0][0]+1) + '\t' + \
                    str(np.where(syntax==st.split('_')[4]+ '_' + st.split('_')[5])[0][0]+1) + '\t0'
    else:
        return '0\t0\t0\t0\t0\t0\t0'

rootDir = "D:\\"
projectDir = os.path.join(rootDir, 'fmri', 'LanguageMVPA')

fnames = os.listdir(os.path.join(projectDir, 'from_scanner'))
subData = defaultdict(list)
[writesubrun(f) for f in fnames]