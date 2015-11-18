#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
# using the chunks to do cross-classification.
from mvpa2.suite import *
import os
import threading
import Queue

# maybe i should detrend...
# INITIALIZING
print "Initializing..."
# initialize paths
projectDir="Z:\\fmri\\LanguageMVPA"
codeDir="Z:\GitHub\LanguageMVPA\multivariate\python"
betaPath = os.path.join(projectDir, "betas")
maskPath = os.path.join(projectDir, "masks", "sub")
labelPath = os.path.join(codeDir, "labels")
outPath = os.path.join(projectDir, "Maps")

# initialize subjects, masks, contrast
contrasts = ["verb", "syntax", "anim", "stimtype", "ActPass", "RelCan", "cross_anim", "cross_verb"]
subList = ["LMVPA001", "LMVPA002", "LMVPA003", "LMVPA005", "LMVPA006", "LMVPA007", "LMVPA008", "LMVPA009", "LMVPA010",
           "LMVPA011", "LMVPA013", "LMVPA014", "LMVPA015", "LMVPA016", "LMVPA017", "LMVPA018", "LMVPA019"]
maskList = ["left_IFG_operc", "IFG_triang", "STG_post", "grayMatter"]
mask = maskList[3]
con = contrasts[2]
if 'cross' in con:
    slType = "cross classification"
    slInt = 0
else:
    slType = "cross validation"
    slInt = 1
# initialize classifier
clf = LinearNuSVMC()
# clf=RbfNuSVMC()
cv = CrossValidation(clf, NFoldPartitioner(), errorfx=lambda p, t: np.mean(p == t),  enable_ca=['stats'])
# 9mm radius searchlight
# sl = sphere_searchlight(cv, radius=4, postproc=mean_sample())
cvSL = sphere_searchlight(cv, radius=4)
attr = SampleAttributes(os.path.join(labelPath, (con + "_attribute_labels.txt")))

# make sure everything is okay.
configMessage = "Searchlight type: " + slType \
    + "\nNumber of subjects: " + str(len(subList))  \
    + "\nMask: " + mask \
    + "\nContrast: " + con

print configMessage

cont = raw_input("Continue? (y/n) \n")
if cont != 'y':
    print "not confirmed, exiting... "
    sys.exit()
else:
    print "Running searchlight analysis... "

# def run_cv_sl(threadName, sl, ds):

def run_cv_sl(sl, ds, t):
    res = sl(ds)
    sphere_accs = res.samples[0]
    map2nifti(ds, sphere_accs).to_filename(os.path.join(outPath, sub + '_' + mask + '_' + con + '_' + t + '_cvsl.nii.gz'))

def run_cc_sl(sl, ds):
    res = sl(ds)
    l2p = res.samples[0]
    p2l = res.samples[1]
    map2nifti(ds, l2p).to_filename(os.path.join(outPath, sub + '_' + mask + '_' + con + '_L2P_ccsl.nii.gz'))
    map2nifti(ds, p2l).to_filename(os.path.join(outPath, sub + '_' + mask + '_' + con + '_P2L_ccsl.nii.gz'))

print "Loaded class descriptions... Let's go!"
for i in range(0, len(subList)):
    sub = subList[i]
    print sub
    bSeriesName = str(sub + "_LSA_Series.nii.gz")
    bSeries = os.path.join(betaPath, bSeriesName)
    maskName = str(sub+"_"+mask+".nii.gz")
    maskFile = os.path.join(maskPath, maskName)
    allFds = fmri_dataset(samples=bSeries, targets=attr.targets, chunks=attr.chunks, mask=maskFile)
    del(bSeries, bSeriesName, maskName, maskFile)

    if slInt == 1:
        lFds = allFds[0:32, :]
        zscore(lFds)
        run_cv_sl(cvSL, lFds, 'lang')
        pFds = allFds[32:64, :]
        zscore(pFds)
        run_cv_sl(cvSL, pFds, 'pic')

    elif slInt == 0:
        zscore(allFds)
        run_cc_sl(cvSL, allFds)
    else:
        print "Something went wrong... exiting. \n"
        sys.exit()


print "Finished!\n"

"""
class CVThread(threading.Thread):
    def __init__(self, threadID, name, sl, ds):
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.name = name
        self.sl = sl
        self.ds = ds

    def run(self):
        print "Starting " + self.name
        run_cv_sl(self.name, self.sl, self.ds)
        print "Exiting " + self.name

thread1 = CVThread(1, "Thread-1", cvSL, lFds)
thread2 = CVThread(2, "Thread-2", cvSL, pFds)
thread1.start()
thread2.start()
"""