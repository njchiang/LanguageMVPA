#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
# this script works the opposite from a traditional pymvpa workflow. it takes pymvpa data,
# pulls the necessary information, runs an encoding model.
import sys
sys.path.append('D:\\GitHub\\LanguageMVPA\\multivariate\\python\\utils')
sys.path.append('/Users/njchiang/GitHub/LanguageMVPA/multivariate/python/utils')
import os
import numpy as np
from mvpa2.misc.io import SampleAttributes
from mvpa2.mappers.zscore import zscore
# import mvpa2.generators.partition as part
# from mvpa2.measures.base import CrossValidation
# from mvpa2.measures.searchlight import sphere_searchlight
# from mvpa2.clfs.skl import SKLLearnerAdapter
# import time
# from mvpa2.misc import errorfx
# local utils
import lmvpautils as lmvpa
from sklearn import linear_model
from mvpa2.mappers.skl_adaptor import SKLTransformer
import ridge as ridge
print "Initializing..."
# initialize paths
paths, subList, contrasts, maskList = lmvpa.initpaths()

# initialize subjects, masks, contrast
mask = "left_STG_post"
con = "syntax"
dsType = "Lang"
# dsType = "Pic"


# load the data
def loadsubdata(m, c, dst):
    if 'cross' in c:
        slt = "cross classification"
        sli = 0
        dst = "Full"
    elif con == "stimtype":
        dst = "Full"
        sli = 1
        slt = "cross validation"
    else:
        slt = "cross validation"
        sli = 1

    # make sure everything is okay.
    configmessage = "Analysis type: " + slt \
        + "\nNumber of subjects: " + str(len(subList))  \
        + "\nMask: " + m \
        + "\nContrast: " + c \
        + "\nSet: " + dst
    print configmessage

    cont = raw_input("Continue? (y/n) \n")
    if cont != 'y':
        print "not confirmed, exiting... "
        sys.exit()

    print "Running " + dst + " dataset"

    cv_attr = SampleAttributes(os.path.join(paths[4], (c + "_attribute_labels.txt")))
    d = []
    for i in range(0, len(subList)):
        sub = subList[i]
        # print sub
        tmp = lmvpa.loadsub(paths, sub, m=mask, a=cv_attr)
        tmp.targets= tmp.targets.astype(float)
        if tmp.shape[1] > 0:
            print "added"
            if dst == "Lang":
                d.append(tmp[0:32])
            elif dst == "Pic":
                d.append(tmp[32:64])
            else:
                d.append(tmp)
    return d


# to run encoder on all voxels, run as a searchlight
def runencoder(d, f):
    # alphas = np.logspace(0, 3, 20)
    alphas=.1
    # run ridgeCV to find best alpha. then apply that alpha and run train on the training set, then test on testing set
    r = []
    for i, n in enumerate(np.unique(d.chunks)):
        c = linear_model.RidgeCV()
        c.fit(f[d.chunks != n], d.samples[d.chunks != n])
        # r.append(c.score(f[d.chunks == n], d.samples[d.chunks == n]))
        r.append(ridge.ridge_corr(
                f[d.chunks != n], f[d.chunks == n], d.samples[d.chunks != n], d.samples[d.chunks == n], c.alpha_))
    return np.array(r)


# load the data
ds_all = loadsubdata(mask, con, dsType)

# try different features-- this does the entire brain region at once... really I want voxel by voxel...
features = np.loadtxt(os.path.join(paths[4], "classFeatures.csv"), skiprows=1, delimiter=',')

_ = [zscore(ds) for ds in ds_all]
res = [runencoder(ds, features[0:len(ds)]) for ds in ds_all]
wsc_results = np.vstack(res)
print np.mean(wsc_results, axis=0)


# # this represent throwing pymvpa2 at ridge regression. This isn't compatible with multidimensional features...
# # initialize the classifier
# # def initCV(nf):
# from mvpa2.clfs.ridge import RidgeReg
# import mvpa2.generators.partition as part
# from mvpa2.measures.base import CrossValidation
# import time
# from mvpa2.misc import errorfx
# def initdecoder():
#     # initialize classifier
#     # this part is deprecated...
#     # rather than predicting voxels based on features.
#     c = RidgeReg()
#     cv = CrossValidation(c, part.NFoldPartitioner(attr='chunks'), errorfx=errorfx.corr_error_prob)
#     return cv
#
# # run within subjects using pymvpa
# def runWSRoi(ds):
#     # inject the subject ID into all datasets
#     # zscore(ds)
#     wsc_start_time = time.time()
#     cv = initdecoder()
#     w = cv(ds)
#     print "done in " + str((time.time() - wsc_start_time,)) + " seconds"
#     # stats for fold
#     return w