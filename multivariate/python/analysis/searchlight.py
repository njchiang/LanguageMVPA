#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
# this script represents just throwing pymvpa at the problem. doesn't work great, and I suspect it's
# because we're using an encoding model.
"""Encoding analysis: using only motion corrected and coregistered data. This analysis applies a Savitsky-Golay filter,
then...
todo: add ridge regression, maybe try crossmodal classification"""
# initialize stuff
import sys
# initialize stuff
if sys.platform == 'darwin':
    plat = 'usb'
    # plat = 'mac'
    sys.path.append('/Users/njchiang/GitHub/LanguageMVPA/multivariate/python/analysis')
    sys.path.append('/Users/njchiang/GitHub/python-fmri-utils/utils')
    debug = True
else:
    plat = 'win'
    sys.path.append('D:\\GitHub\\LanguageMVPA\\multivariate\\python\\analysis')
    sys.path.append('D:\\GitHub\\python-fmri-utils\\utils')
    debug = False
import lmvpautils as lmvpa
from sklearn import linear_model as lm
from mvpa2.clfs.skl import SKLLearnerAdapter
from mvpa2.clfs import svm
# debug = True
thisContrast = 'verb'
roi = 'grayMatter'
filterLen = 49
filterOrd = 3
r = 4 #searchlight radius
clf = svm.LinearNuSVMC()
paths, subList, contrasts, maskList = lmvpa.initpaths(plat)
if debug:
    subList = {'LMVPA005': subList['LMVPA005']}
mc_params = lmvpa.loadmotionparams(paths, subList)
beta_events = lmvpa.loadevents(paths, subList)

import numpy as np
import os
from mvpa2.datasets.mri import map2nifti
from mvpa2.mappers.zscore import zscore
import SavGolFilter as sg
import mvpa2.datasets.eventrelated as er
from mvpa2.measures.base import CrossValidation
from mvpa2.generators.partition import NFoldPartitioner
import searchlightutils as sl

for sub in subList.keys():
    thisSub = {sub: subList[sub]}
    dsdict = lmvpa.loadsubdata(paths, thisSub, m=roi, c='trial_type')
    thisDS = dsdict[sub]
    mc_params = lmvpa.loadmotionparams(paths, thisSub)
    beta_events = lmvpa.loadevents(paths, thisSub)
    # savitsky golay filtering
    sg.sg_filter(thisDS, filterLen, filterOrd)
    # gallant group zscores before regression.

    # zscore w.r.t. rest trials
    # zscore(thisDS, param_est=('targets', ['rest']), chunks_attr='chunks')
    # zscore entire set. if done chunk-wise, there is no double-dipping (since we leave a chunk out at a time).
    zscore(thisDS, chunks_attr='chunks')
    print "beta extraction"
    ## BETA EXTRACTION ##
    rds, events = lmvpa.amendtimings(thisDS.copy(), beta_events[sub])
    evds = er.fit_event_hrf_model(rds, events, time_attr='time_coords',
                                  condition_attr=('trial_type', 'chunks'),
                                  design_kwargs={'add_regs': mc_params[sub], 'hrf_model': 'canonical'},
                                  return_model=True)

    fds = lmvpa.replacetargets(evds, contrasts, thisContrast)
    fds = fds[fds.targets != '0']

    print "searchlights"
    ## initialize classifier
    cv = CrossValidation(clf, NFoldPartitioner())
    from mvpa2.measures.searchlight import sphere_searchlight
    cvSL = sphere_searchlight(cv, radius=r)


    # now I have betas per chunk. could just correlate the betas, or correlate the predictions for corresponding runs
    lidx = fds.chunks < fds.sa['chunks'].unique[len(fds.sa['chunks'].unique)/2]
    pidx = fds.chunks >= fds.sa['chunks'].unique[len(fds.sa['chunks'].unique) / 2]

    lres = sl.run_cv_sl(cvSL, fds[lidx].copy(deep=False))
    pres = sl.run_cv_sl(cvSL, fds[pidx].copy(deep=False))
    from mvpa2.base import dataset
    map2nifti(thisDS, dataset.vstack([lres, pres])).\
        to_filename(os.path.join(
                    paths[0], 'Maps', 'PyMVPA',
                    sub + '_' + roi + '_' + thisContrast + '_cvsl.nii.gz'))

    del lres, pres, cvSL

    cvSL = sphere_searchlight(cv, radius=r)
    crossSet = fds.copy()
    crossSet.chunks[lidx] = 1
    crossSet.chunks[pidx] = 2
    cres = sl.run_cv_sl(cvSL, crossSet.copy(deep=False))
    map2nifti(thisDS, cres[0]).to_filename(
        os.path.join(paths[0], 'Maps', 'PyMVPA',
                     sub + '_' + roi + '_' + (thisContrast) + '_P2L.nii.gz'))
    map2nifti(thisDS, cres[1]).to_filename(
        os.path.join(paths[0], 'Maps', 'PyMVPA',
                     sub + '_' + roi + '_' + (thisContrast) + '_L2P.nii.gz'))

