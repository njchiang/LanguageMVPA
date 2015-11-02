#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
#configure for windows

from mvpa2.suite import *
import os

projectDir="Z:\\fmri\\LanguageMVPA"
mask="grayMatter"
constrasts=["verb", "syntax", "stimtype"]
con="stimtype"

subList = ["LMVPA001", "LMVPA002", "LMVPA003", "LMVPA005", "LMVPA006", "LMVPA007", "LMVPA008", "LMVPA009", "LMVPA010", "LMVPA011", "LMVPA013", "LMVPA014", "LMVPA015", "LMVPA016", "LMVPA017", "LMVPA018", "LMVPA019"]
betaPath=os.path.join(projectDir, "betas")
maskPath=os.path.join(projectDir, "masks", "sub")
labelPath=os.path.join(projectDir, "labels")
maskList=["IFG_operc", "IFG_triang", "STG_post"]
mask="grayMatter"
clf=LinearNuSVMC()
#clf=RbfNuSVMC()
fsel = SensitivityBasedFeatureSelection(OneWayAnova(), FixedNElementTailSelector(500, mode='select', tail='upper'))
fclf = FeatureSelectionClassifier(clf, fsel)
#cv = CrossValidation(fclf, NFoldPartitioner(),errorfx=lambda p, t: np.mean(p == t),  enable_ca=['stats'])
cv = CrossValidation(clf, NFoldPartitioner(),errorfx=lambda p, t: np.mean(p == t),  enable_ca=['stats'])

unNormAcc=np.zeros(len(subList))
normAcc=np.zeros(len(subList))
attr=SampleAttributes(os.path.join(labelPath, (con+"_attribute_labels.txt")))	

for i in range(0, len(subList)):
	sub=subList[i]
	print sub
	bSeriesName=str(sub + "_LSA_Series.nii.gz")
	maskName=str(sub+"_"+mask+".nii.gz")

	bSeries=os.path.join(betaPath, bSeriesName)
	maskFile=os.path.join(maskPath, maskName)
	
	allFds=fmri_dataset(samples=bSeries, targets=attr.targets, chunks=attr.chunks, mask=maskFile)

	oResults=cv(allFds)
	print np.round(cv.ca.stats.stats['ACC%'],1)
	print cv.ca.stats.matrix
	unNormAcc[i]=np.round(cv.ca.stats.stats['ACC%'],1)

	zFds=allFds
	zscore(zFds, chunks_attr='chunks')
	zResults=cv(zFds)
	print np.round(cv.ca.stats.stats['ACC%'],1)
	print cv.ca.stats.matrix

	normAcc[i]=np.round(cv.ca.stats.stats['ACC%'],1)
	

print np.mean(unNormAcc)
print np.mean(normAcc)
"""
	oResults=cv(allFds)
	print np.round(cv.ca.stats.stats['ACC%'],1)
	print cv.ca.stats.matrix

	zFds=allFds
	zscore(zFds, chunks_attr='chunks')
	zResults=cv(zFds)
	print np.round(cv.ca.stats.stats['ACC%'],1)
	

	lFds=allFds[0:32,:]
	lResults=cv(lFds)
	print np.round(cv.ca.stats.stats['ACC%'],1)
	print cv.ca.stats.matrix
	unNormAcc[i]=np.round(cv.ca.stats.stats['ACC%'],1)
	
	zscore(lFds, chunks_attr=None)
	lzResults=cv(lFds)
	print np.round(cv.ca.stats.stats['ACC%'],1)	
	print cv.ca.stats.matrix
	normAcc[i]=np.round(cv.ca.stats.stats['ACC%'],1)
	
	pFds=allFds[32:64,:]
	pResults=cv(pFdsFds)
	print np.round(cv.ca.stats.stats['ACC%'],1)
	print cv.ca.stats.matrix
	unNormAcc[i]=np.round(cv.ca.stats.stats['ACC%'],1)
	
	zscore(pFds, chunks_attr=None)
	pzResults=cv(pFds)
	print np.round(cv.ca.stats.stats['ACC%'],1)	
	print cv.ca.stats.matrix
	normAcc[i]=np.round(cv.ca.stats.stats['ACC%'],1)
"""

"""The final preprocessing step is data-normalization. This is a required step
for many classification algorithms. It scales all features (voxels)
into approximately the same range and removes the mean. In this example, we
perform a chunk-wise normalization and compute standard deviation and mean for
z-scoring based on the volumes corresponding to rest periods in the experiment.
The resulting features could be interpreted as being voxel salience relative
to 'rest'.

zscore(dataset, chunks_attr='chunks', param_est=('targets', ['rest']), dtype='float32')

After normalization is completed, we no longer need the 'rest'-samples and
remove them."""

#dataset = dataset[dataset.sa.targets != 'rest']

