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
    plat = 'usb'
    # plat = 'mac'
    sys.path.append('/Users/njchiang/GitHub/LanguageMVPA/multivariate/python/utils')
    debug = True
else:
    plat = 'win'
    sys.path.append('D:\\GitHub\\LanguageMVPA\\multivariate\\python\\utils')
    debug = False
import lmvpautils as lmvpa
debug = False
thisContrast = ['syntax', 'verb']
roi = 'grayMatter'
filterLen = 49
filterOrd = 3
alpha=150
paths, subList, contrasts, maskList = lmvpa.initpaths(plat)
if debug:
    subList = {'LMVPA005': subList['LMVPA005']}
# load things in as trial type for easy regression, then swap out labels accordingly
# do we actually want to load all data simultaneously? or one brain at a time...
# ds_all = lmvpa.loadsubdata(paths, subList, m=roi, c='trial_type')
# motion parameters for all subjects
mc_params = lmvpa.loadmotionparams(paths, subList)
# events for beta extraction
# beta_events = lmvpa.loadevents(paths, subList, c='trial_type')
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


# i really need to just train on chunks though...
# add optimization...
def encodingcorr(betas, ds, idx=None, part_attr='chunks'):
    # iterate through the attributes of the dataset
    # get the betas that correspond to this.. but how if we've already picked?
    if not idx is None:
        des = np.array(betas.sa['regressors']).T[idx,:]
        betas = betas
        ds = ds[idx].copy()
    else:
        des = np.array(betas.sa['regressors']).T
    # need to only pull the correct betas...
    res = []
    for i in ds.sa[part_attr].unique:
        trainidx = ds.sa['chunks'].unique[ds.sa['chunks'].unique != i]
        thesebetas = []
        for j in trainidx:
            thesebetas.append(betas.samples[betas.chunks==j])
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

