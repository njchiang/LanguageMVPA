#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
# configure for windows
# whoa. rerunning same CV ruins it for future iterations. weird.
from mvpa2.suite import *
import os
print "Initializing..."
# initialize paths
projectDir="Z:\\fmri\\LanguageMVPA"
codeDir="Z:\GitHub\LanguageMVPA\multivariate\python"
# betaPath = os.path.join(projectDir, "betas", "tstat")
# betaPath = os.path.join(projectDir, "betas", "orig")
betaPath = os.path.join(projectDir, "betas", "tstat")
maskPath = os.path.join(projectDir, "masks", "sub")
labelPath = os.path.join(codeDir, "labels")
outPath = os.path.join(projectDir, "Maps")

# initialize subjects, masks, contrast
contrasts = ["verb", "syntax", "anim", "stimtype", "ActPass", "RelCan", "cross_anim", "cross_verb"]
# subList = ["LMVPA001", "LMVPA002", "LMVPA003", "LMVPA005", "LMVPA006", "LMVPA007", "LMVPA008", "LMVPA009", "LMVPA010",
#            "LMVPA011", "LMVPA013", "LMVPA014", "LMVPA015", "LMVPA016", "LMVPA017", "LMVPA018", "LMVPA019"]
# subList = ["LMVPA001", "LMVPA002", "LMVPA003", "LMVPA005", "LMVPA006", "LMVPA008", "LMVPA009", "LMVPA010",
#           "LMVPA011", "LMVPA013", "LMVPA014", "LMVPA015", "LMVPA016", "LMVPA017", "LMVPA018", "LMVPA019"]
subList = ["testv2"]
maskList = ["left_IFG_operc", "left_IFG_triang", "left_STG_post", "left_MTG_post", "grayMatter"]
mask = maskList[4]
con = contrasts[2]
if 'cross' in con:
    slType = "cross classification"
    slInt = 0
else:
    slType = "cross validation"
    slInt = 1
# initialize classifier
# clf = LinearNuSVMC()
clf = LinearCSVMC()
# clf=RbfNuSVMC()
nVox = 100
fsel = SensitivityBasedFeatureSelection(OneWayAnova(), FixedNElementTailSelector(nVox, mode='select', tail='upper'))
fclf = FeatureSelectionClassifier(clf, fsel)
cv = CrossValidation(fclf, NFoldPartitioner(), errorfx=lambda p, t: np.mean(p == t),  enable_ca=['stats'])
# cv = CrossValidation(fclf, OddEvenPartitioner(), errorfx=lambda p, t: np.mean(p == t),  enable_ca=['stats'])
# cv = CrossValidation(clf, NFoldPartitioner(), errorfx=lambda p, t: np.mean(p == t),  enable_ca=['stats'])
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
