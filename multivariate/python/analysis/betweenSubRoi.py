#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
# configure for windows
from mvpa2.suite import *
import os
print "Initializing..."
# initialize paths
projectDir = "D:\\fmri\\LanguageMVPA"
codeDir = "D:\GitHub\LanguageMVPA\multivariate\python"
betaType = "crossSub"  # tstat cope orig
betaPath = os.path.join(projectDir, "betas", betaType)
maskPath = os.path.join(projectDir, "masks")
labelPath = os.path.join(codeDir, "labels")
outPath = os.path.join(projectDir, "Maps")

# initialize subjects, masks, contrast
contrasts = ["verb", "syntax", "anim"]
maskList = ["left_IFG_operc", "left_IFG_triang", "left_STG_post", "left_MTG_post", "langNet", "grayMatter"]
mask = maskList[5]
con = contrasts[2]


def loadData(stim, c, m):
    cv_attr = SampleAttributes(os.path.join(labelPath, ("crossSub_" + c + "_attribute_labels.txt")))
    bSeriesName = str(stim + "_" + c + ".nii.gz")
    maskName = str("2mm_" + m + ".nii.gz")
    bSeries = os.path.join(betaPath, bSeriesName)
    maskFile = os.path.join(maskPath, maskName)
    ds = (fmri_dataset(samples=bSeries, targets=cv_attr.targets, chunks=cv_attr.chunks, mask=maskFile))
    return ds

lds = loadData("lang", con, mask)
pds = loadData("pic", con, mask)

def initCV():
    # initialize classifier
    # clf = LinearCSVMC()
    clf = LinearNuSVMC()
    # clf = RbfNuSVMC()
    # feature selection helpers
    nf = 100
    fselector = FixedNElementTailSelector(nf, tail='upper',
                                          mode='select',sort=False)
    sbfs = SensitivityBasedFeatureSelection(OneWayAnova(), fselector,
                                            enable_ca=['sensitivities'])
    # create classifier with automatic feature selection
    fsclf = FeatureSelectionClassifier(clf, sbfs)
    cv = CrossValidation(fsclf,
                         NFoldPartitioner(attr='chunks'),
                         errorfx=mean_match_accuracy)
    # cv = CrossValidation(fclf, NFoldPartitioner(), errorfx=lambda p, t: np.mean(p == t),  enable_ca=['stats'])
    return cv

lcv = initCV()
bsc_l_results = lcv(lds)
pcv = initCV()
bsc_p_results = pcv(pds)
print "Avg between-subject accuracy: " + str(np.mean(bsc_l_results)) + " +/- " + str(np.std(bsc_l_results) / np.sqrt(16))
print "Avg between-subject accuracy: " + str(np.mean(bsc_p_results)) + " +/- " + str(np.std(bsc_p_results) / np.sqrt(16))

ccv = initCV()
cc = str("cross_" + con)
cv_attr = SampleAttributes(os.path.join(labelPath, ("crossSub_" + cc + "_attribute_labels.txt")))
bSeriesName = str(con + ".nii.gz")
maskName = str("2mm_" + mask + ".nii.gz")
bSeries = os.path.join(betaPath, bSeriesName)
maskFile = os.path.join(maskPath, maskName)
fds = (fmri_dataset(samples=bSeries, targets=cv_attr.targets, chunks=cv_attr.chunks, mask=maskFile))
bsc_c_results = ccv(fds)
print "Avg between-subject accuracy: " + str(np.mean(bsc_c_results[0]))
print "Avg between-subject accuracy: " + str(np.mean(bsc_c_results[1]))
