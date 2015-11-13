#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
#configure for windows
# simple searchlight script...
from mvpa2.suite import *
import os

projectDir="Z:\\fmri\\LanguageMVPA"
codeDir="Z:\GitHub\LanguageMVPA\multivariate\python"
mask="grayMatter"
constrasts=["verb", "syntax", "anim", "stimtype"]
con="anim"

subList = ["LMVPA001", "LMVPA002", "LMVPA003", "LMVPA005", "LMVPA006", "LMVPA007", "LMVPA008", "LMVPA009", "LMVPA010", "LMVPA011", "LMVPA013", "LMVPA014", "LMVPA015", "LMVPA016", "LMVPA017", "LMVPA018", "LMVPA019"]
betaPath=os.path.join(projectDir, "betas")
maskPath=os.path.join(projectDir, "masks", "sub")
labelPath=os.path.join(projectDir, "labels")
outPath=os.path.join(projectDir, "Maps")
maskList=["IFG_operc", "IFG_triang", "STG_post"]
mask="grayMatter"
clf=LinearNuSVMC()
#clf=RbfNuSVMC()
cv = CrossValidation(clf, NFoldPartitioner(),errorfx=lambda p, t: np.mean(p == t),  enable_ca=['stats'])
#9mm radius searchlight
sl = sphere_searchlight(cv, radius=3, postproc=mean_sample())

attr=SampleAttributes(os.path.join(labelPath, (con+"_attribute_labels.txt")))	


for i in range(0, len(subList)):
    sub=subList[i]
    print sub
    bSeriesName=str(sub + "_LSA_Series.nii.gz")
    maskName=str(sub+"_"+mask+".nii.gz")
    bSeries=os.path.join(betaPath, bSeriesName)
    maskFile=os.path.join(maskPath, maskName)
    allFds=fmri_dataset(samples=bSeries, targets=attr.targets, chunks=attr.chunks, mask=maskFile)
    lFds=allFds[0:32,:]
    zscore(lFds, chunks_attr=None)
    lRes=sl(lFds)
    sphere_accs = lRes.samples[0]
    map2nifti(lFds, sphere_accs).to_filename(os.path.join(outPath, sub+'_' + con + 'sl.nii.gz'))


