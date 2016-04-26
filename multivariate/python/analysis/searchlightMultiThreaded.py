#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
# using the chunks to do cross-classification.
import sys
sys.path.append('D:\\GitHub\\LanguageMVPA\\multivariate\\python\\analysis')
from mvpa2.suite import *
import lmvpautils as lmvpa
import os

paths, subList, contrasts, maskList = lmvpa.initpaths()

# initialize subjects, masks, contrast
# mask = maskList[4]
# con = contrasts[5]
mask = "grayMatter"
con = "RelCan"
if 'cross' in con:
    slType = "cross classification"
    slInt = 0
else:
    slType = "cross validation"
    slInt = 1

# initialize labels and classifier
# if len(con) > 1:
#     print "Running RSA"
#     clf = lmvpa.initrsa(con)
# else:
#     print "Running SVM"
attr = SampleAttributes(os.path.join(paths[4], (con + "_attribute_literal_labels.txt")))
clf = LinearNuSVMC()
# clf=RbfNuSVMC()
cv = CrossValidation(clf, NFoldPartitioner())

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


def error2acc(d):
    d.samples *= -1
    d.samples += 1
    return d


def run_cv_sl(sl, ds, t, p):
    fds = ds.copy(deep=False, sa=['targets', 'chunks'], fa=['voxel_indices'], a=['mapper'])
    # is this necessary if i'm doing tstats? in toolbox it says zscore w.r.t. rest condition,
    # but technically a tstat is already normalized?
    zscore(fds)
    wsc_start_time = time.time()
    print "running " + str(t) + " at " + time.strftime("%H:%M:%S")
    thisSL = sl
    res = thisSL(fds)
    print "done in " + str((time.time() - wsc_start_time,)) + " seconds"
    res = error2acc(res)
    map2nifti(res, imghdr=ds.a.imghdr).to_filename(os.path.join(p[0], 'Maps', 'PyMVPA', sub + '_' + mask + '_' + con + '_' + t + '_cvsl.nii.gz'))


def run_cc_sl(sl, ds, p):
    fds = ds.copy(deep=False, sa=['targets', 'chunks'], fa=['voxel_indices'], a=['mapper'])
    zscore(fds)
    thisSL = sl
    res = thisSL(fds)
    res = error2acc(res)
    l2p = res.samples[0]
    p2l = res.samples[1]
    map2nifti(ds, l2p).to_filename(os.path.join(p[0], 'Maps', 'PyMVPA', sub + '_' + mask + '_' + con + '_L2P_ccsl.nii.gz'))
    map2nifti(ds, p2l).to_filename(os.path.join(p[0], 'Maps', 'PyMVPA',, sub + '_' + mask + '_' + con + '_P2L_ccsl.nii.gz'))

print "Loaded class descriptions... Let's go!"
for i in range(0, len(subList)):
    sub = subList[i]
    fullSet = lmvpa.loadsubbetas(paths, sub, m=mask, a=attr)
    if slInt == 1:
        cvSL = sphere_searchlight(cv, radius=2, postproc=mean_sample())
        lFds = fullSet[0:32, :]
        run_cv_sl(cvSL, lFds, 'lang', paths)
        # run_splitcv_sl(cvSL, lFds, 'lang')
        pFds = fullSet[32:64, :]
        run_cv_sl(cvSL, pFds, 'pic', paths)
        # run_splitcv_sl(cvSL, pFds, 'pic')

    elif slInt == 0:
        cvSL = sphere_searchlight(cv, radius=2)
        run_cc_sl(cvSL, fullSet, paths)
    else:
        print "Something went wrong... exiting. \n"
        sys.exit()


print "Finished!\n"

