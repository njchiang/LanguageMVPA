#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
# configure for windows
from mvpa2.suite import *
import os
import threading

projectDir = "Z:\\fmri\\LanguageMVPA"
codeDir = "Z:\GitHub\LanguageMVPA\multivariate\python"
constrasts = ["verb", "syntax", "stimtype", "anim"]
con = "syntax"

subList = ["LMVPA001", "LMVPA002", "LMVPA003", "LMVPA005", "LMVPA006", "LMVPA007", "LMVPA008", "LMVPA009",
           "LMVPA010", "LMVPA011", "LMVPA013", "LMVPA014", "LMVPA015", "LMVPA016", "LMVPA017", "LMVPA018", "LMVPA019"]
# subList = ["LMVPA017"]
betaPath = os.path.join(projectDir, "betas")
maskPath = os.path.join(projectDir, "masks", "sub")
labelPath = os.path.join(projectDir, "labels")
maskList = ["IFG_operc", "IFG_triang", "STG_post"]
mask = "left_IFG_operc"
clf = LinearNuSVMC()
# clf = RbfNuSVMC()
# fsel = SensitivityBasedFeatureSelection(OneWayAnova(), FixedNElementTailSelector(500, mode='select', tail='upper'))
# fclf = FeatureSelectionClassifier(clf, fsel)
# cv = CrossValidation(fclf, NFoldPartitioner(),errorfx=lambda p, t: np.mean(p == t),  enable_ca=['stats'])
cv = CrossValidation(clf, NFoldPartitioner(), errorfx=lambda p, t: np.mean(p == t), enable_ca=['stats'])

langAcc = np.zeros(len(subList))
picAcc = np.zeros(len(subList))
cv_attr = SampleAttributes(os.path.join(labelPath, (con + "_attribute_labels.txt")))
# cross_attr = SampleAttributes(os.path.join(labelPath, ("cross_"+con+"_attribute_labels.txt")))
for i in range(0, len(subList)):
    sub = subList[i]
    print sub
    bSeriesName = str(sub + "_LSA_Series.nii.gz")
    maskName = str(sub + "_" + mask + ".nii.gz")
    bSeries = os.path.join(betaPath, bSeriesName)
    maskFile = os.path.join(maskPath, maskName)
    print "loading files..."
    allFds = fmri_dataset(samples=bSeries, targets=cv_attr.targets, chunks=cv_attr.chunks, mask=maskFile)
    lFds = allFds[0:32,:]
    # inanimLFds=lFds[0:16,:]
    animLFds=lFds[16:32,:]
    lResults = cv(animLFds)
    print np.round(cv.ca.stats.stats['ACC%'], 1)
    print cv.ca.stats.matrix
    langAcc[i] = np.round(cv.ca.stats.stats['ACC%'], 1)
    pFds = allFds[32:64,:]
    # inanimLFds=lFds[0:16,:]
    animPFds=lFds[16:32,:]
    pResults = cv(animPFds)
    print np.round(cv.ca.stats.stats['ACC%'], 1)
    print cv.ca.stats.matrix
    picAcc[i] = np.round(cv.ca.stats.stats['ACC%'], 1)
    # zscore(zFds, chunks_attr='chunks')
print "Language Mean: " + str(np.mean(langAcc)) + " Standard deviation: " + str(np.std(langAcc))
print "Picture Mean: " + str(np.mean(picAcc)) + " Standard deviation: " + str(np.std(picAcc))


