#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
# this script represents just throwing pymvpa at the problem. doesn't work great, and I suspect it's
# because we're using an encoding model.
"""
Current version of encoding analysis: takes timing files for betas and makes design matrix out of them.
Two ways to do cross-validation:
    1) run same regressors on each chunk
    2) run everything at once, but have a separate set of regressors per chunk and index them accordingly.
This script runs version (1)
let R = number of regressors
betas will be 4R by nVox.
"""
import sys
# initialize stuff
if sys.platform == 'darwin':
    plat = 'usb'
    # plat = 'mac'
    sys.path.append('/Users/njchiang/GitHub/python-fmri-utils/utils')
    debug = True
else:
    plat = 'win'
    sys.path.append('D:\\GitHub\\python-fmri-utils\\utils')
    debug = False

import lmvpautils as lmvpa
debug = False
thisContrast = ['syntax']
roi = 'grayMatter'
filterLen = 49
filterOrd = 3
alpha = 1000

paths, subList, contrasts, maskList = lmvpa.initpaths(plat)
if debug:
    subList = {'LMVPA005': subList['LMVPA005']}

# ds_all = lmvpa.loadsubdata(paths, subList, m=roi, c='trial_type')
# motion parameters for all subjects
mc_params = lmvpa.loadmotionparams(paths, subList)
# events for beta extraction
# add everything as a sample attribute
beta_events = lmvpa.loadevents(paths, subList)

import sklearn.linear_model as lm
import SKLMapper as sklm
import BootstrapRidgeMapper as bsr
import numpy as np
import os
from mvpa2.datasets.mri import map2nifti
from mvpa2.mappers.zscore import zscore
import SavGolFilter as sg

# use the betas to predict timeseries for left out chunk
def encodingcorr(betas, ds, idx=None, part_attr='chunks'):
    # iterate through the attributes of the dataset
    # get the betas that correspond to this.. but how if we've already picked?
    if not idx is None:
        des = np.array(betas.sa['regressors']).T[idx, :]
        ds = ds[idx].copy()
    else:
        des = np.array(betas.sa['regressors']).T
    # need to only pull the correct betas...
    res = []
    for i in ds.sa[part_attr].unique:
        trainidx = ds.sa['chunks'].unique[ds.sa['chunks'].unique != i]
        thesebetas = []
        for j in trainidx:
            thesebetas.append(betas.samples[betas.chunks == j])
        # estbetas = np.vstack(thesebetas)
        estbetas = np.mean(np.dstack(thesebetas), axis=-1) # take mean of all betas to predict...
        pred = np.dot(des[np.array(ds.sa[part_attr]) == i][:, np.array(betas.sa[part_attr]) == i],
                      estbetas)


        # model_params.sa[part] = np.repeat(ds.sa[part].unique[c],
        #                                   len(model_params), axis=0)
        # mds.append(model_params)
        resvar = (ds.samples[ds.chunks == i] - pred).var(0)
        Rsqs = 1 - (resvar / ds.samples[ds.chunks == i].var(0))
        corrs = np.sqrt(np.abs(Rsqs)) * np.sign(Rsqs)
        res.append(corrs)
    from mvpa2.datasets import Dataset
    return Dataset(np.vstack(res), sa={part_attr: ds.sa[part_attr].unique})


def findalpha(ds, X, part='chunks', nchunks=2, nboots=50, alphas=None):
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

    # we can model out motion and just not use those betas.
    # Ridge
    if isinstance(thisContrast, basestring):
        thisContrast=[thisContrast]
    # do i want to do this crossed with chunks for everything?...
    desX, rds = lmvpa.make_designmat(rds, events, time_attr='time_coords', condition_attr=thisContrast,
                                     design_kwargs={'hrf_model': 'canonical', 'drift_model': 'blank'},
                                     regr_attrs=None)

    searchstring = 'glm_label' # only change if you want to do a specific subset
    regressor_names = []
    for rn in rds.sa.keys():
        if searchstring in rn:
            regressor_names.append(rn)
    regressor_names.sort()
    regressor_names.append('constant')
    # alpha = np.logspace(0, 3, 20)
    # chunks refers to the sa. seems to be a copying method.
    # language within
    lidx = thisDS.chunks < thisDS.sa['chunks'].unique[len(thisDS.sa['chunks'].unique)/2]
    lclf = sklm.SKLRegressionMapper(regs=regressor_names, add_regs=[], clf=lm.Ridge(alpha=alpha), part_attr='chunks',
                                    return_design=True)
    betas = lclf(rds) # not sure if i need to regress out chunkwise mean too
    lbetas = betas.copy()
    lres = encodingcorr(betas, thisDS, lidx, part_attr='chunks')
    # now I have betas per chunk. could just correlate the betas, or correlate the predictions for corresponding runs
    print 'language ' + str(np.mean(lres))
    # map2nifti(thisDS, np.mean(lres, axis=0)).to_filename(
    #     os.path.join(paths[0], 'Maps', 'Encoding', sub + '_' + roi + '_' + '+'.join(thisContrast) + '_Lridge.nii.gz'))

    # pictures within
    pidx = thisDS.chunks >= thisDS.sa['chunks'].unique[len(thisDS.sa['chunks'].unique) / 2]
    pres = encodingcorr(betas, thisDS, pidx, part_attr='chunks')
    print 'pictures: ' + str(np.mean(pres))
    # map2nifti(thisDS, np.mean(pres, axis=0)).to_filename(
    #     os.path.join(paths[0], 'Maps', 'Encoding', sub + '_' + roi + '_' + '+'.join(thisContrast) + '_Pridge.nii.gz'))
    from mvpa2.base import dataset
    map2nifti(thisDS, dataset.vstack([lres, pres])) \
        .to_filename(os.path.join(paths[0], 'Maps', 'Encoding', sub + '_' + roi + '_' + '+'.join(thisContrast) +
                                  '_a' + str(alpha) +'_ridge.nii.gz'))
    del lres,pres
    crossSet = thisDS.copy()
    crossSet.chunks[lidx] = 1
    crossSet.chunks[pidx] = 2
    cres = encodingcorr(betas, crossSet, part_attr='chunks')
    print 'cross: ' + str(np.mean(cres))
    map2nifti(thisDS, cres[0]).to_filename(
        os.path.join(paths[0], 'Maps', 'Encoding', sub + '_' + roi + '_' + '+'.join(thisContrast) + '_P2Lridge.nii.gz'))
    map2nifti(thisDS, cres[1]).to_filename(
        os.path.join(paths[0], 'Maps', 'Encoding',sub + '_' + roi + '_' + '+'.join(thisContrast) + '_L2Pridge.nii.gz'))

