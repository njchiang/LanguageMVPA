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
projectDir="D:\\fmri\\LanguageMVPA"
codeDir="D:\GitHub\LanguageMVPA\multivariate\python"
betaPath = os.path.join(projectDir, "betas", "cope")
maskPath = os.path.join(projectDir, "masks", "sub")
labelPath = os.path.join(codeDir, "labels")
outPath = os.path.join(projectDir, "Maps", "PyMVPA")

# initialize subjects, masks, contrast
contrasts = ["verb", "syntax", "anim", "stimtype", "ActPass", "RelCan", "cross_anim", "cross_verb"]
# subList = ["LMVPA001", "LMVPA002", "LMVPA003", "LMVPA005", "LMVPA006", "LMVPA007", "LMVPA008", "LMVPA009", "LMVPA010",
#            "LMVPA011", "LMVPA013", "LMVPA014", "LMVPA015", "LMVPA016", "LMVPA017", "LMVPA018", "LMVPA019"]
subList = ["LMVPA019", "LMVPA018", "LMVPA017", "LMVPA016", "LMVPA015", "LMVPA014", "LMVPA013", "LMVPA011", "LMVPA010",
           "LMVPA009", "LMVPA008", "LMVPA007", "LMVPA006", "LMVPA005", "LMVPA003", "LMVPA002", "LMVPA001"]
# subList = ["LMVPA001"]
maskList = ["left_IFG_operc", "left_IFG_triang", "left_STG_post", "left_MTG_post", "grayMatter"]
mask = maskList[4]
con = contrasts[4]
if 'cross' in con:
    slType = "cross classification"
    slInt = 0
else:
    slType = "cross validation"
    slInt = 1
# initialize classifier
clf = LinearNuSVMC()
# clf=RbfNuSVMC()
# cv = CrossValidation(clf, NFoldPartitioner(), errorfx=lambda p, t: np.mean(p == t),  enable_ca=['stats'])
# cv = CrossValidation(clf, NFoldPartitioner(), errorfx=lambda p, t: np.mean(p == t))
# note: this calculates error, and flips is back.
cv = CrossValidation(clf, NFoldPartitioner())
# 9mm radius searchlight
# sl = sphere_searchlight(cv, radius=4, postproc=mean_sample())
# cvSL = sphere_searchlight(cv, radius=4)
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

# need to be fixed, I think it writes all of the samples instead of the mean


def error2acc(d):
    d.samples *= -1
    d.samples += 1
    return d

def run_cv_sl(sl, ds, t):
    fds = ds.copy(deep=False, sa=['targets', 'chunks'], fa=['voxel_indices'], a=['mapper'])
    zscore(ds)
    print "running " + str(t)
    thisSL = sl
    res = thisSL(ds)
    print "done with " + str(t)
    res = error2acc(res)
    map2nifti(res, imghdr=ds.a.imghdr).to_filename(os.path.join(outPath, sub + '_' + mask + '_' + con + '_' + t + '_cvsl.nii.gz'))

# these need to be fixed... I think it writes all of the samples instead of the mean
def run_splitcv_sl(sl, ds, t):
    fds = ds.copy(deep=False, sa=['targets', 'chunks'], fa=['voxel_indices'], a=['mapper'])
    ids = fds[0:16, :]
    zscore(ids)
    thisSL = sl
    res = thisSL(ids)
    res = error2acc(res)
    map2nifti(res, imghdr=ds.a.imghdr).to_filename(os.path.join(outPath, sub + '_' + mask + '_inanim_' + con + '_' + t + '_cvsl.nii.gz'))
    ads = ds[16:32, :]
    zscore(ads)
    thisSL = sl
    res = thisSL(ads)
    res = error2acc(res)
    map2nifti(res, imghdr=ds.a.imghdr).to_filename(os.path.join(outPath, sub + '_' + mask + '_anim_' + con + '_' + t + '_cvsl.nii.gz'))

def run_cc_sl(sl, ds):
    fds = ds.copy(deep=False, sa=['targets', 'chunks'], fa=['voxel_indices'], a=['mapper'])
    zscore(fds)
    thisSL = sl
    res = thisSL(fds)
    res = error2acc(res)
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
    fullSet = fmri_dataset(samples=bSeries, targets=attr.targets, chunks=attr.chunks, mask=maskFile)
    del(bSeries, bSeriesName, maskName, maskFile)

    if slInt == 1:
        cvSL = sphere_searchlight(cv, radius=2, postproc=mean_sample())
        lFds = fullSet[0:32, :]
        run_cv_sl(cvSL, lFds, 'lang')
        # run_splitcv_sl(cvSL, lFds, 'lang')
        pFds = fullSet[32:64, :]
        run_cv_sl(cvSL, pFds, 'pic')
        # run_splitcv_sl(cvSL, pFds, 'pic')

    elif slInt == 0:
        cvSL = sphere_searchlight(cv, radius=2)
        run_cc_sl(cvSL, fullSet)
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