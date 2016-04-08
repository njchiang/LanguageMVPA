#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
import sys
sys.path.append('D:\\GitHub\\LanguageMVPA\\multivariate\\python\\utils')
sys.path.append('/Users/njchiang/GitHub/LanguageMVPA/multivariate/python/utils')
import mvpa2.clfs as clfs
from mvpa2.measures.base import CrossValidation
from mvpa2.misc.io import SampleAttributes
import lmvpautils as lmvpa
import ridge as ridge
import os
import numpy as np
from sklearn import linear_model

print "Initializing..."
# initialize paths
paths, subList, contrasts, maskList = lmvpa.initpaths()

# nVox = 100
# initialize subjects, masks, contrast
mask = "left_STG_post"
con = "cross_anim"
# dsType = "Lang"
dsType = "Pic"


# initialize the classifier
# def initCV(nf):
def initdecoder():
    # initialize classifier
    # this part is deprecated...
    # rather than predicting voxels based on features.
    c = clfs.ridge.RidgeReg()
    cv = CrossValidation(c, clfs.meta.NFoldPartitioner(attr='chunks'), errorfx=corr_error_prob)
    return cv

def runencoder(d,f):
    alphas = np.logspace(0, 3, 20)
    ridge.ridge_corr(d.targets, d.targes, d.samples, d.samples, alphas)
    # determine the best alpha...
    # c = linear_model.Ridge() sklearn option

# load the data
def loadSubData(m, c, t):
    cv_attr = SampleAttributes(os.path.join(paths[4], (c + "_attribute_labels.txt")))
    d = []
    for i in range(0, len(subList)):
        sub = subList[i]
        # print sub
        tmp = lmvpa.loadsub(paths, sub, m=mask, a=cv_attr)
        tmp.targets= tmp.targets.astype(float)
        if tmp.shape[1] > 0:
            print "added"
            if dsType == "Lang":
                d.append(tmp[0:32])
            elif dsType == "Pic":
                d.append(tmp[32:64])
            else:
                d.append(tmp)
    return d


# run the analysis
def runWSRoi(ds):
    # inject the subject ID into all datasets
    # zscore(ds)
    wsc_start_time = time.time()
    cv = initdecoder()
    w = cv(ds)
    print "done in " + str((time.time() - wsc_start_time,)) + " seconds"
    # stats for fold
    return w

if 'cross' in con:
    slType = "cross classification"
    slInt = 0
    dsType = "Full"
elif con == "stimtype":
    dsType = "Full"
    sInt = 1
    slType = "cross validation"
else:
    slType = "cross validation"
    slInt = 1

# make sure everything is okay.
configMessage = "Analysis type: " + slType \
    + "\nNumber of subjects: " + str(len(subList))  \
    + "\nMask: " + mask \
    + "\nContrast: " + con \
    + "\nSet: " + dsType
print configMessage

cont = raw_input("Continue? (y/n) \n")
if cont != 'y':
    print "not confirmed, exiting... "
    sys.exit()

# load the data
ds_all = loadSubData(mask, con, dsType)

print "Running " + dsType + " dataset"
res = [runWSRoi(d) for d in ds_all]
wsc_results = np.vstack(res)

"""plan
make a new ds_copy that has the voxel activity as "targets", and the features as samples. then can train...
"""