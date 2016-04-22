#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
# this script represents just throwing pymvpa at the problem. doesn't work great, and I suspect it's
# because we're using an encoding model.
"""Encoding analysis: using only motion corrected and coregistered data. This analysis applies a Savitsky-Golay filter,
then...
todo: add ridge regression, maybe try crossmodal classification"""
import sys
# initialize stuff
if sys.platform == 'darwin':
    plat = 'mac'
    sys.path.append('/Users/njchiang/GitHub/LanguageMVPA/multivariate/python/utils')
    debug = True
else:
    plat = 'win'
    sys.path.append('D:\\GitHub\\LanguageMVPA\\multivariate\\python\\utils')
    debug = False
import lmvpautils as lmvpa
# plat = 'usb'
# debug = True
paths, subList, contrasts, maskList = lmvpa.initpaths(plat)
thisContrast = 'syntax'
roi = 'langNet'
filterLen = 49
filterOrd = 3
if debug:
    subList = {'LMVPA003': subList['LMVPA003']}
# load things in as trial type for easy regression, then swap out labels accordingly
# do we actually want to load all data simultaneously? or one brain at a time...
# ds_all = lmvpa.loadsubdata(paths, subList, m=roi, c='trial_type')
# motion parameters for all subjects
mc_params = lmvpa.loadmotionparams(paths, subList)
# events for beta extraction
# beta_events = lmvpa.loadevents(paths, subList, c='trial_type')
# add everything as a sample attribute
beta_events = lmvpa.loadevents(paths, subList)

# import sklearn.linear_model as lm
# import SKLMapper as sklm
import BootstrapRidgeMapper as bsr
import numpy as np
import copy
import os
from mvpa2.datasets.mri import map2nifti
from mvpa2.mappers.zscore import zscore
import SavGolFilter as sg
import mvpa2.datasets.eventrelated as er

# i really need to just train on chunks though...
def runCVBootstrap(ds, X, part='chunks', nchunks=2, nboots=100, alphas=None):
    # runs cross validation on the chunks of the dataset (leave-one-out)
    if alphas is None:
        alphas = np.logspace(0, 3, 20)
    # push design into source dataset
    # mds=[]
    res = []
    from mvpa2.datasets import Dataset
    for c in np.arange(len(ds.sa[part].unique)):
        # need to combine regressors into one...
        trainidx = ds.sa[part].value != ds.sa[part].unique[c]
        testidx = ds.sa[part].value == ds.sa[part].unique[c]
        glm_regs = [(reg, X.matrix[trainidx, i]) for i, reg in enumerate(X.names)]
        # MUST CONTAIN  chunklen and nchunks
        glmfit_kwargs = {'chunklen': np.sum(trainidx) / 5 / nchunks,
                         'nchunks': nchunks,
                         'alphas': alphas,
                         'nboots': nboots,
                         'corrmin': 0.2,
                         'joined': None,
                         'singcutoff': 1e-10,
                         'normalpha': False,
                         'single_alpha': True,
                         'use_corr': True}
        # first parameter is null because we have made our own design matrix (add_regs)
        # now this is working. need to add functionality to append chunks
        bootstrapridge = bsr.BootstrapRidgeMapper([], glmfit_kwargs=glmfit_kwargs,
                                                  add_regs=glm_regs,
                                                  return_model=True)

        model_params = bootstrapridge(ds[trainidx])
        pred = np.dot(X.matrix[testidx], model_params.samples)
        # model_params.sa[part] = np.repeat(ds.sa[part].unique[c],
        #                                   len(model_params), axis=0)
        # mds.append(model_params)
        resvar = (ds.samples[testidx] - pred).var(0)
        Rsqs = 1 - (resvar / ds.samples[testidx].var(0))
        corrs = np.sqrt(np.abs(Rsqs)) * np.sign(Rsqs)
        res.append(corrs)
    return Dataset(np.vstack(res), sa={part: ds.sa[part].unique})
    # later this will loop


for sub in subList.keys():
    thisSub = {sub: subList[sub]}
    dsdict = lmvpa.loadsubdata(paths, thisSub, m=roi, c='trial_type')
    thisDS = dsdict[sub]
    # savitsky golay filtering
    sg.sg_filter(thisDS, filterLen, filterOrd)
    # gallant group zscores before regression.

    # zscore w.r.t. rest trials
    # zscore(thisDS, param_est=('targets', ['rest']), chunks_attr='chunks')
    # zscore entire set. if done chunk-wise, there is no double-dipping (since we leave a chunk out at a time).
    zscore(thisDS, chunks_attr='chunks')

    # kay method: leave out a model run, use it to fit an HRF for each voxel
    # huth method: essentially use FIR
    # mumford method: deconvolution with canonical HRF

    # refit events and regress...
    # get timing data from timing files
    rds, events = lmvpa.amendtimings(thisDS.copy(), beta_events[sub])
    # estimate HRF from picture runs for language, and vice versa
    # now: make timing files for each feature that encompass every timepoint but with different intensities
    # desX = lmvpa.make_designmat(rds, events, time_attr='time_coords', condition_attr=['targets', 'chunks'],
    #                             design_kwargs={'add_regs': mc_params[sub], 'hrf_model': 'canonical'},
    #                             glmfit_kwargs=None, regr_attrs=None)
    # GLM
    # normal regression. doesn't use desX from above.
    evds = er.fit_event_hrf_model(rds, events, time_attr='time_coords',
                                  condition_attr=(thisContrast, 'chunks'),
                                  design_kwargs={'add_regs': mc_params[sub], 'hrf_model': 'canonical'},
                                  return_model=True)

    # Ridge
    desX, rds = lmvpa.make_fulldesignmat(rds, events, time_attr='time_coords', condition_attr=[thisContrast],
                                     design_kwargs={'hrf_model': 'canonical', 'drift_model': 'blank'},
                                     glmfit_kwargs=None, regr_attrs=None)
    # language within
    lidx = thisDS.chunks < thisDS.sa['chunks'].unique[len(thisDS.sa['chunks'].unique)/2]
    lds = copy.copy(desX)
    lds.matrix = lds.matrix[lidx]
    lres = runCVBootstrap(rds.copy()[lidx], lds)
    print 'language ' + str(np.mean(lres))
    map2nifti(thisDS, np.mean(lres, axis=0)).to_filename(os.path.join(paths[0], 'Maps', 'Encoding', sub + '_' + roi + '_' + thisContrast + '_Lridge.nii.gz'))
    del lds, lres  # just cleaning up
    # pictures within
    pidx = thisDS.chunks >= thisDS.sa['chunks'].unique[len(thisDS.sa['chunks'].unique) / 2]
    pds = copy.copy(desX)
    pds.matrix = pds.matrix[pidx]
    pres = runCVBootstrap(rds.copy()[pidx], pds)
    print 'pictures: ' + str(np.mean(pres))
    map2nifti(thisDS, np.mean(pres, axis=0)).to_filename(os.path.join(paths[0], 'Maps', 'Encoding', sub + '_' + roi + '_' + thisContrast + '_Pridge.nii.gz'))
    del pds, pres

    crossSet = rds.copy()
    crossSet.chunks[lidx] = 1
    crossSet.chunks[pidx] = 2
    cres = runCVBootstrap(crossSet, desX)
    print 'cross: ' + str(np.mean(cres))
    map2nifti(thisDS, cres[0]).to_filename(os.path.join(paths[0], 'Maps', 'Encoding', sub + '_' + roi + '_' + thisContrast + '_P2Lridge.nii.gz'))
    map2nifti(thisDS, cres[1]).to_filename(os.path.join(paths[0], 'Maps', 'Encoding', sub + '_' + roi + '_' + thisContrast + '_L2Pridge.nii.gz'))



