#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
# configure for windows
from mvpa2.suite import *
import os
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
maskList = ["left_IFG_operc", "left_IFG_triang", "left_STG_post", "left_MTG_post", "grayMatter"]
mask = maskList[0]
con = contrasts[1]
if 'cross' in con:
    slType = "cross classification"
    slInt = 0
else:
    slType = "cross validation"
    slInt = 1
# initialize classifier
clf = LinearNuSVMC()
# clf=RbfNuSVMC()

fsel = SensitivityBasedFeatureSelection(OneWayAnova(), FixedNElementTailSelector(500, mode='select', tail='upper'))
fclf = FeatureSelectionClassifier(clf, fsel)
cv = CrossValidation(fclf, NFoldPartitioner(), errorfx=lambda p, t: np.mean(p == t),  enable_ca=['stats'])
#cv = CrossValidation(clf, NFoldPartitioner(), errorfx=lambda p, t: np.mean(p == t),  enable_ca=['stats'])
attr = SampleAttributes(os.path.join(labelPath, (con + "_attribute_labels.txt")))

# make sure everything is okay.
configMessage = "Analysis type: " + slType \
    + "\nNumber of subjects: " + str(len(subList))  \
    + "\nMask: " + mask \
    + "\nContrast: " + con

print configMessage

cont = raw_input("Continue? (y/n) \n")
if cont != 'y':
    print "not confirmed, exiting... "
    sys.exit()
else:
    print "Running ROI analysis... "

langAcc = np.zeros(len(subList))
picAcc = np.zeros(len(subList))
animLangAcc = np.zeros(len(subList))
animPicAcc = np.zeros(len(subList))
inAnimLangAcc = np.zeros(len(subList))
inAnimPicAcc = np.zeros(len(subList))
cv_attr = SampleAttributes(os.path.join(labelPath, (con + "_attribute_labels.txt")))
for i in range(0, len(subList)):
    sub = subList[i]
    print sub
    bSeriesName = str(sub + "_LSA_Series.nii.gz")
    maskName = str(sub + "_" + mask + ".nii.gz")
    bSeries = os.path.join(betaPath, bSeriesName)
    maskFile = os.path.join(maskPath, maskName)
    print "loading files..."
    allFds = fmri_dataset(samples=bSeries, targets=cv_attr.targets, chunks=cv_attr.chunks, mask=maskFile)
    lFds = allFds[0:32, :]

    print "running language trials: inanimate"
    inAnimLFds = lFds[0:16, :]
    zscore(inAnimLFds)
    ilResults = cv(inAnimLFds)
    del(inAnimLFds)
    print np.round(cv.ca.stats.stats['ACC%'], 1)
    print cv.ca.stats.matrix
    inAnimLangAcc[i] = np.round(cv.ca.stats.stats['ACC%'], 1)

    print "running language trials: animate"
    animLFds = lFds[16:32, :]
    zscore(animLFds)
    alResults = cv(animLFds)
    del(animLFds)
    print np.round(cv.ca.stats.stats['ACC%'], 1)
    print cv.ca.stats.matrix
    animLangAcc[i] = np.round(cv.ca.stats.stats['ACC%'], 1)

    print " running all language trials"
    zscore(lFds)
    lResults = cv(lFds)
    print np.round(cv.ca.stats.stats['ACC%'], 1)
    print cv.ca.stats.matrix
    langAcc[i] = np.round(cv.ca.stats.stats['ACC%'], 1)


    pFds = allFds[32:64, :]
    print "running picture trials: inanimate"
    inAnimPFds = pFds[0:16, :]
    zscore(inAnimPFds)
    ipResults = cv(inAnimPFds)
    del(inAnimPFds)
    print np.round(cv.ca.stats.stats['ACC%'], 1)
    print cv.ca.stats.matrix
    inAnimPicAcc[i] = np.round(cv.ca.stats.stats['ACC%'], 1)

    print "running picture trials: animate"
    animPFds = pFds[16:32, :]
    zscore(animPFds)
    apResults = cv(animPFds)
    del(animPFds)
    print np.round(cv.ca.stats.stats['ACC%'], 1)
    print cv.ca.stats.matrix
    animPicAcc[i] = np.round(cv.ca.stats.stats['ACC%'], 1)

    print "running all picture trials"
    zscore(pFds)
    pResults = cv(pFds)
    print np.round(cv.ca.stats.stats['ACC%'], 1)
    print cv.ca.stats.matrix
    picAcc[i] = np.round(cv.ca.stats.stats['ACC%'], 1)


print "Inanim Language Mean: " + str(np.mean(inAnimLangAcc)) + " Standard deviation: " + str(np.std(inAnimLangAcc))
print "Anim Language Mean: " + str(np.mean(animLangAcc)) + " Standard deviation: " + str(np.std(animLangAcc))
print "Language Mean: " + str(np.mean(langAcc)) + " Standard deviation: " + str(np.std(langAcc))
print "Inanim Picture Mean: " + str(np.mean(inAnimPicAcc)) + " Standard deviation: " + str(np.std(inAnimPicAcc))
print "Anim Picture Mean: " + str(np.mean(animPicAcc)) + " Standard deviation: " + str(np.std(animPicAcc))
print "Picture Mean: " + str(np.mean(picAcc)) + " Standard deviation: " + str(np.std(picAcc))


# can split by language/picture (duh), then split by anim/inanim if i need...
