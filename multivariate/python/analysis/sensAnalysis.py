#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

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
subList = ["LMVPA001", "LMVPA002", "LMVPA003", "LMVPA005", "LMVPA006", "LMVPA008", "LMVPA009", "LMVPA010",
           "LMVPA011", "LMVPA013", "LMVPA014", "LMVPA015", "LMVPA016", "LMVPA017", "LMVPA018", "LMVPA019"]
# subList = ["testv2"]
maskList = ["left_IFG_operc", "left_IFG_triang", "left_STG_post", "left_MTG_post", "grayMatter"]
mask = maskList[0]
con = contrasts[1]
sub = "LMVPA001"
print sub
bSeriesName = str(sub + "_LSA_Series.nii.gz")
maskName = str(sub + "_" + mask + ".nii.gz")
bSeries = os.path.join(betaPath, bSeriesName)
maskFile = os.path.join(maskPath, maskName)
print "loading files..."
cv_attr = SampleAttributes(os.path.join(labelPath, (con + "_attribute_labels.txt")))

allFds = fmri_dataset(samples=bSeries, targets=cv_attr.targets, chunks=cv_attr.chunks, mask=maskFile)
sensanas = {
    'a) ANOVA': OneWayAnova(postproc=absolute_features()),
    'b) Linear SVM weights': LinearNuSVMC().get_sensitivity_analyzer(
                                               postproc=absolute_features()),
    'c) Splitting ANOVA (odd-even)':
        RepeatedMeasure(
            OneWayAnova(postproc=absolute_features()),
            OddEvenPartitioner()),
    'd) Splitting SVM (odd-even)':
        RepeatedMeasure(
            LinearNuSVMC().get_sensitivity_analyzer(postproc=absolute_features()),
            OddEvenPartitioner()),
    'e) Splitting ANOVA (nfold)':
        RepeatedMeasure(
            OneWayAnova(postproc=absolute_features()),
            NFoldPartitioner()),
    'f) Splitting SVM (nfold)':
        RepeatedMeasure(
            LinearNuSVMC().get_sensitivity_analyzer(postproc=absolute_features()),
            NFoldPartitioner()),
            }
"""    'h) GNB Searchlight':
        sphere_gnbsearchlight(GNB(), NFoldPartitioner(cvtype=1),
                              radius=0, errorfx=mean_match_accuracy)
     'c) I-RELIEF': IterativeRelief(postproc=absolute_features()),
     """



fig = 0
pl.figure(figsize=(14, 8))

keys = sensanas.keys()
keys.sort()

for s in keys:
    # tell which one we are doing
    print "Running %s ..." % (s)

    sana = sensanas[s]
    # compute sensitivies
    sens = sana(allFds)
    # I-RELIEF assigns zeros, which corrupts voxel masking for pylab's
    # imshow, so adding some epsilon :)
    smap = sens.samples[0] + 0.00001

    # map sensitivity map into original dataspace
    orig_smap = allFds.mapper.reverse1(smap)
    masked_orig_smap = np.ma.masked_array(orig_smap, mask=orig_smap == 0)

    # make a new subplot for each classifier
    fig += 1
    pl.subplot(3, 3, fig)

    pl.title(s)

    pl.imshow(masked_orig_smap[..., 0].T,
             interpolation='nearest',
             aspect=1.25,
             cmap=pl.cm.autumn)

    # uniform scaling per base sensitivity analyzer
    ## if s.count('ANOVA'):
    ##     pl.clim(0, 30)
    ## elif s.count('SVM'):
    ##     pl.clim(0, 0.055)
    ## else:
    ##     pass

    pl.colorbar(shrink=0.6)