#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
import sys
sys.path.append('D:\\GitHub\\LanguageMVPA\\multivariate\\python\\analysis')
from mvpa2.suite import *
import lmvpautils as lmvpa
import os

print "Initializing..."
# initialize paths
paths, subList, contrasts, maskList = lmvpa.initpaths()
subList = ["LMVPA001", "LMVPA002", "LMVPA005", "LMVPA006", "LMVPA008", "LMVPA009", "LMVPA010",
           "LMVPA011", "LMVPA013", "LMVPA014", "LMVPA015", "LMVPA016", "LMVPA017", "LMVPA018", "LMVPA019"]
# nVox = 100
# initialize subjects, masks, contrast
mask = "lang_semantics_left_IFG"
con = "verb"
dsType = "Lang"
# dsType = "Pic"


# initialize the classifier
# def initCV(nf):
def initCV():
    # initialize classifier
    # clf = LinearCSVMC()
    clf = LinearNuSVMC()
    # clf = RbfNuSVMC()
    # feature selection helpers*
    # fselector = FixedNElementTailSelector(nf, tail='upper',
    #                                       mode='select',sort=False)
    # sbfs = SensitivityBasedFeatureSelection(OneWayAnova(), fselector,
    #                                         enable_ca=['sensitivities'])
    # create classifier with automatic feature selection
    # fsclf = FeatureSelectionClassifier(clf, sbfs)
    # cv = CrossValidation(fsclf,
    cv = CrossValidation(clf,
                         NFoldPartitioner(attr='chunks'),
                         errorfx=mean_match_accuracy)
    # cv = CrossValidation(fclf, NFoldPartitioner(), errorfx=lambda p, t: np.mean(p == t),  enable_ca=['stats'])
    return cv


# load the data
def loadSubData(m, c, t):
    # subFileName = os.path.join(betaPath, t + "_" + m + "_" + c + ".nii.gz")
    cv_attr = SampleAttributes(os.path.join(paths[4], (c + "_attribute_labels.txt")))
    d = []
    for i in range(0, len(subList)):
        sub = subList[i]
        # print sub
        tmp = lmvpa.loadsub(paths, sub, m=mask, a=cv_attr)
        if dsType == "Lang":
            d.append(tmp[0:32])
        elif dsType == "Pic":
            d.append(tmp[32:64])
        else:
            d.append(tmp)
    # h5save(subFileName, d)
    return d


# run the analysis
def runWSRoi(fullDataset):
    # inject the subject ID into all datasets
    for i, sd in enumerate(fullDataset):
        sd.sa['subject'] = np.repeat(i, len(sd))
    _ = [zscore(ds) for ds in fullDataset]
    wsc_start_time = time.time()
    cv = initCV()
    wsc_results = [cv(j) for j in fullDataset]
    wsc_results = vstack(wsc_results)
    print "done in " + str((time.time() - wsc_start_time,)) + " seconds"
    print "Avg within-subject accuracy: " + str(np.mean(wsc_results)) + " +/- " + str(np.std(wsc_results) / np.sqrt(len(fullDataset) - 1))
    # stats for fold
    for f in range(0, len(fullDataset[0].UC)):
        tmp = np.take(wsc_results.samples, np.arange(f, len(wsc_results), len(fullDataset[0].UC)))
        print "Avg for fold " + str(f+1) + ": " + str(np.mean(tmp)) + " +/- " + str(np.std(tmp) / np.sqrt(len(tmp) - 1))
        del tmp
    # stats for sub
    for s in range(0, len(fullDataset)):
        tmp = np.take(wsc_results.samples, np.arange(s*len(fullDataset[0].UC), (s+1)*len(fullDataset[0].UC)))
        print "Avg for sub " + str(s+1) + ": " + str(np.mean(tmp)) + " +/- " + str(np.std(tmp) / np.sqrt(len(tmp) - 1))
        del tmp
    return wsc_results

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

# number of subjects
print "\n\nDataset Review..."
print "Number of subjects: " + str(len(ds_all))
# number of categories
print "Number of categories: " + str(len(ds_all[0].UT))
# number of runs
print "Number of folds: " + str(len(ds_all[0].UC))
# number timepoints
print "Number of timepoints: " + str(len(ds_all[0]))
# numVoxels
print "Number of voxels: " + str((ds_all[0].shape[1]))
cont = raw_input("Continue? (y/n) \n")
if cont != 'y':
    print "not confirmed, exiting... "
    sys.exit()
else:
    print "Running ROI analysis... "

print "Running " + dsType + " dataset"
res = runWSRoi(ds_all)


"""
  try:
        print "Found previously generated file, loading..."
        d = h5load(subFileName)
    except IOError:
        print "Could not load dataset, regenerating..." """
