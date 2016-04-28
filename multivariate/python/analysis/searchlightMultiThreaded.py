#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
# using the chunks to do cross-classification.
import sys
sys.path.append('D:\\GitHub\\LanguageMVPA\\multivariate\\python\\utils')
from mvpa2.suite import *
import lmvpautils as lmvpa
import os
import searchlightutils as sl
paths, subList, contrasts, maskList = lmvpa.initpaths('win')

# initialize subjects, masks, contrast
# mask = maskList[4]
# con = contrasts[5]
mask = "grayMatter"
con = "syntax"
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
attr = SampleAttributes(os.path.join(paths[1], "labels", (con + "_attribute_literal_labels.txt")))
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

print "Loaded class descriptions... Let's go!"
for i in range(0, len(subList)):
    sub = subList.keys()[i]
    fullSet = lmvpa.loadsubbetas(paths, sub, m=mask, a=attr)
    if slInt == 1:
        cvSL = sphere_searchlight(cv, radius=2)
        lFds = fullSet[0:32, :]
        lres = sl.run_cv_sl(cvSL, lFds, 'lang', paths)
        # run_splitcv_sl(cvSL, lFds, 'lang')
        pFds = fullSet[32:64, :]
        pres = sl.run_cv_sl(cvSL, pFds, 'pic', paths)
        # run_splitcv_sl(cvSL, pFds, 'pic')
        map2nifti(fullSet, vstack([lres, pres])).\
            to_filename(os.path.join(
            paths[0], 'Maps', 'PyMVPA', sub + '_' + mask + '_' + con + '_cvsl.nii.gz'))

    elif slInt == 0:
        cvSL = sphere_searchlight(cv, radius=2)
        l2p, p2l = sl.run_cc_sl(cvSL, fullSet, paths)
        # ERROR HERE, THEY'RE BACKWARDS
        map2nifti(fullSet, l2p).to_filename(
            os.path.join(
                paths[0], 'Maps', 'PyMVPA', sub + '_' + mask + '_' + con + '_L2P_ccsl.nii.gz'))
        map2nifti(fullSet, p2l).to_filename(
            os.path.join(
                paths[0], 'Maps', 'PyMVPA', sub + '_' + mask + '_' + con + '_P2L_ccsl.nii.gz'))

    else:
        print "Something went wrong... exiting. \n"
        sys.exit()


print "Finished!\n"

