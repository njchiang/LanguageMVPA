#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
import numpy as np
import pylab as pl
from os.path import join as pjoin
from mvpa2.suite import *
import os
from os.path import join as pjoin

print "Initializing..."

projectDir="D:\\fmri\\LanguageMVPA"
codeDir="D:\GitHub\LanguageMVPA\multivariate\python"
betaType = "tstat"  # tstat cope orig
betaPath = os.path.join(projectDir, "betas", betaType)
maskPath = os.path.join(projectDir, "masks", "sub")
labelPath = os.path.join(codeDir, "labels")
outPath = os.path.join(projectDir, "Maps")

# initialize subjects, masks, contrast
contrasts = ["verb", "syntax", "anim",  "ActPass", "RelCan", "stimtype", "cross_anim", "cross_verb"]
subList = ["LMVPA001", "LMVPA002", "LMVPA003", "LMVPA005", "LMVPA006", "LMVPA007", "LMVPA008", "LMVPA009", "LMVPA010",
           "LMVPA011", "LMVPA013", "LMVPA014", "LMVPA015", "LMVPA016", "LMVPA017", "LMVPA018", "LMVPA019"]
# subList = ["LMVPA001", "LMVPA002", "LMVPA003", "LMVPA005", "LMVPA006", "LMVPA008", "LMVPA009", "LMVPA010",
#            "LMVPA011", "LMVPA013", "LMVPA014", "LMVPA015", "LMVPA016", "LMVPA017", "LMVPA018", "LMVPA019"]
subList = ["LMVPA005"]
maskList = ["left_IFG_operc", "left_IFG_triang", "left_STG_post", "left_MTG_post", "langNet", "grayMatter"]
mask = maskList[0]
con = contrasts[1]
dsType = "Lang"
# dsType = "Pic"

def loadSubData(m, c, t):
    # subFileName = os.path.join(betaPath, t + "_" + m + "_" + c + ".nii.gz")
    cv_attr = SampleAttributes(os.path.join(labelPath, (c + "_attribute_labels.txt")))
    d = []
    for i in range(0, len(subList)):
        sub = subList[i]
        print sub
        bSeriesName = str(sub + "_LSA_Series.nii.gz")
        maskName = str(sub + "_" + mask + ".nii.gz")
        bSeries = os.path.join(betaPath, bSeriesName)
        maskFile = os.path.join(maskPath, maskName)
        print "loading files..."
        return fmri_dataset(samples=bSeries, targets=cv_attr.targets, chunks=cv_attr.chunks, mask=maskFile)
        # tmp = (fmri_dataset(samples=bSeries, targets=cv_attr.targets, chunks=cv_attr.chunks, mask=maskFile))
    #     if dsType == "Lang":
    #         d.append(tmp[0:32])
    #     elif dsType == "Pic":
    #         d.append(tmp[32:64])
    #     else:
    #         d.append(tmp)
    # # h5save(subFileName, d)
    # return d

def plot_mtx(mtx, labels, title):
    pl.figure()
    pl.imshow(mtx, interpolation='nearest')
    pl.xticks(range(len(mtx)), labels, rotation=-45)
    pl.yticks(range(len(mtx)), labels)
    pl.title(title)
    pl.clim((0,1))
    pl.colorbar()


ds_all = loadSubData(mask, con, dsType)
ds = ds_all[0:32]
mtgs = mean_group_sample(['targets'])
mtds = mtgs(ds)

from mvpa2.measures import rsa
dsm = rsa.PDist(square=True)
res = dsm(mtds)
plot_mtx(res, mtds.sa.targets, 'ROI pattern correlation distances')

from mvpa2.measures.searchlight import sphere_searchlight
dsm = rsa.PDist(square=False)
sl = sphere_searchlight(dsm, 2)
slres = sl(mtds)

# score each searchlight sphere result wrt global pattern dissimilarity
distinctiveness = np.sum(np.abs(slres), axis=0)
print 'Most dissimilar patterns around', \
        mtds.fa.voxel_indices[distinctiveness.argmax()]
# take a look at the this dissimilarity structure
from scipy.spatial.distance import squareform
plot_mtx(squareform(slres.samples[:, distinctiveness.argmax()]),
         mtds.sa.targets,
         'Maximum distinctive searchlight pattern correlation distances')

mtcgs = mean_group_sample(['targets', 'chunks'])
mtcds = mtcgs(ds)

# searchlight consistency measure -- how correlated are the structures
# across runs
dscm = rsa.PDistConsistency()
sl_cons = sphere_searchlight(dscm, 2)
slres_cons = sl_cons(mtcds)

# mean correlation
mean_consistency = np.mean(slres_cons, axis=0)
print 'Most stable dissimilarity patterns around', \
        mtds.fa.voxel_indices[mean_consistency.argmax()]
# Look at this pattern
plot_mtx(squareform(slres.samples[:, mean_consistency.argmax()]),
         mtds.sa.targets,
         'Most consistent searchlight pattern correlation distances')


""" THIS IS THE ONE WE WANT """
# let's see where in the brain we find dissimilarity structures that are
# similar to our most stable one
tdsm = rsa.PDistTargetSimilarity(
            slres.samples[:, mean_consistency.argmax()])
# using a searchlight
from mvpa2.base.learner import ChainLearner
from mvpa2.mappers.shape import TransposeMapper
sl_tdsm = sphere_searchlight(ChainLearner([tdsm, TransposeMapper()]), 2)
slres_tdsm = sl_tdsm(mtds)


# plot the spatial distribution using NiPy
vol = ds.a.mapper.reverse1(slres_tdsm.samples[0])
import nibabel as nb
anat = nb.load(os.path.join(betaPath, str("LMVPA005" + "_LSA_Series.nii.gz")))

from nipy.labs.viz_tools.activation_maps import plot_map
pl.figure(figsize=(15,4))
sp = pl.subplot(121)
pl.title('Distribution of target similarity structure correlation')
slices = plot_map(
            vol,
            ds.a.imgaffine,
            cut_coords=np.array((12,-42,-20)),
            threshold=.5,
            cmap="bwr",
            vmin=0,
            vmax=1.,
            axes=sp,
            anat=anat.get_data(),
            anat_affine=anat.get_affine(),
         )
img = pl.gca().get_images()[1]
cax = pl.axes([.05, .05, .05, .9])
pl.colorbar(img, cax=cax)

sp = pl.subplot(122)
pl.hist(slres_tdsm.samples[0],
        #range=(0,410),
        normed=False,
        bins=30,
        color='0.6')
pl.ylabel("Number of voxels")
pl.xlabel("Target similarity structure correlation")