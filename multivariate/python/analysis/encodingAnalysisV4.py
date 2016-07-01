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
# bootstrapping works. now need to fix design matrix
betas will be 4R by nVox.
should i include probe in the testing? can zero out all betas for probes
compare to regular regression as baseline
should each probe be individually coded? yes.
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
    sys.path.append('D:\\GitHub\\LanguageMVPA\\multivariate\\python\\analysis')
    sys.path.append('D:\\GitHub\\python-fmri-utils\\utils')
    debug = False

import lmvpautils as lmvpa
debug = False
# thisContrast = ['anim', 'verb', 'ap', 'cr']
thisContrast = ['anim', 'verb', 'ap', 'cr']
roi = 'grayMatter'
filterLen = 49
filterOrd = 3
chunklen=10
nchunks = 2

paths, subList, contrasts, maskList = lmvpa.initpaths(plat)
if debug:
    subList = {'LMVPA005': subList['LMVPA005']}

# ds_all = lmvpa.loadsubdata(paths, subList, m=roi, c='trial_type')
# motion parameters for all subjects
mc_params = lmvpa.loadmotionparams(paths, subList)
# events for beta extraction
# add everything as a sample attribute
beta_events = lmvpa.loadevents(paths, subList)

import BootstrapRidgeMapper as bsr
import numpy as np
import os
from mvpa2.datasets.mri import map2nifti
from mvpa2.mappers.zscore import zscore
import SavGolFilter as sg
import BootstrapRidge as bsr
import copy as cp


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

    # univariate
    # import mvpa2.datasets.eventrelated as er
    #
    # evds = er.fit_event_hrf_model(rds, events, time_attr='time_coords',
    #                               condition_attr=(thisContrast, 'chunks'),
    #                               design_kwargs={'add_regs': mc_params[sub], 'hrf_model': 'canonical'},
    #                               return_model=True)

    # we can model out motion and just not use those betas.
    # Ridge
    if isinstance(thisContrast, basestring):
        thisContrast = [thisContrast]
    # instead of binarizing each one, make them parametric
    desX, rds = lmvpa.make_designmat(rds, events, time_attr='time_coords', condition_attr=thisContrast,
                                     design_kwargs={'hrf_model': 'canonical', 'drift_model': 'blank'},
                                     regr_attrs=None)
    # want to collapse ap and cr, but have anim separate
    des = lmvpa.make_parammat(desX)

    # set chunklen and nchunks
    # split by language and pictures
    lidx = thisDS.chunks < thisDS.sa['chunks'].unique[len(thisDS.sa['chunks'].unique) / 2]
    pidx = thisDS.chunks >= thisDS.sa['chunks'].unique[len(thisDS.sa['chunks'].unique) / 2]
    ldes = cp.copy(des)
    pdes = cp.copy(des)

    ldes.matrix = ldes.matrix[lidx]
    pdes.matrix = pdes.matrix[pidx]
    lwts, lalphas, lres = bsr.bootstrap_ridge(rds[lidx], ldes, chunklen=chunklen, nchunks=nchunks,
                                              part_attr='chunks', mode='test',
                                              alphas=None, single_alpha=True, normalpha=False,
                                              nboots=15, corrmin=.2, singcutoff=1e-10, joined=None,
                                              use_corr=True)
    pwts, palphas, pres = bsr.bootstrap_ridge(rds[pidx], pdes, chunklen=chunklen, nchunks=nchunks,
                                              part_attr='chunks', mode='test',
                                              alphas=None, single_alpha=True, normalpha=False,
                                              nboots=15, corrmin=.2, singcutoff=1e-10, joined=None,
                                              use_corr=True)
    # now I have betas per chunk. could just correlate the betas, or correlate the predictions for corresponding runs
    print 'language ' + str(np.mean(lres))
    # map2nifti(thisDS, np.mean(lres, axis=0)).to_filename(
    #     os.path.join(paths[0], 'Maps', 'Encoding', sub + '_' + roi + '_' + '+'.join(thisContrast) + '_Lridge.nii.gz'))

    # pictures within
    print 'pictures: ' + str(np.mean(pres))
    # map2nifti(thisDS, np.mean(pres, axis=0)).to_filename(
    #     os.path.join(paths[0], 'Maps', 'Encoding', sub + '_' + roi + '_' + '+'.join(thisContrast) + '_Pridge.nii.gz'))
    from mvpa2.base import dataset
    map2nifti(thisDS, dataset.vstack([lres, pres])) \
        .to_filename(os.path.join(paths[0], 'Maps', 'Encoding', sub + '_' + roi + '_' + '+'.join(thisContrast) +
                                  '_cvAlpha_ridge.nii.gz'))
    map2nifti(thisDS, dataset.vstack([lwts, pwts])) \
        .to_filename(os.path.join(paths[0], 'Maps', 'Encoding', sub + '_' + roi + '_' + '+'.join(thisContrast) +
                                  '_cvAlpha_weights.nii.gz'))
    map2nifti(thisDS, dataset.vstack([lres, pres])) \
        .to_filename(os.path.join(paths[0], 'Maps', 'Encoding', sub + '_' + roi + '_' + '+'.join(thisContrast) +
                                  '_cvAlpha_alphas.nii.gz'))

    del lres, pres
    crossSet = thisDS.copy()
    crossSet.chunks[lidx] = 1
    crossSet.chunks[pidx] = 2
    cwts, calphas, cres = bsr.bootstrap_ridge(crossSet, des, chunklen=chunklen, nchunks=2*nchunks,
                                              part_attr='chunks', mode='test',
                                              alphas=None, single_alpha=True, normalpha=False,
                                              nboots=15, corrmin=.2, singcutoff=1e-10, joined=None,
                                              use_corr=True)
    print 'cross: ' + str(np.mean(cres))
    map2nifti(thisDS, cres[0]).to_filename(
        os.path.join(paths[0], 'Maps', 'Encoding', sub + '_' + roi + '_' + '+'.join(thisContrast) + '_P2Lridge.nii.gz'))
    map2nifti(thisDS, cres[1]).to_filename(
        os.path.join(paths[0], 'Maps', 'Encoding', sub + '_' + roi + '_' + '+'.join(thisContrast) + '_L2Pridge.nii.gz'))

    map2nifti(thisDS, cwts[0]).to_filename(
        os.path.join(paths[0], 'Maps', 'Encoding', sub + '_' + roi + '_' + '+'.join(thisContrast) + '_P2Lweights.nii.gz'))
    map2nifti(thisDS, cwts[1]).to_filename(
        os.path.join(paths[0], 'Maps', 'Encoding', sub + '_' + roi + '_' + '+'.join(thisContrast) + '_L2Pweights.nii.gz'))

    map2nifti(thisDS, calphas[0]).to_filename(
        os.path.join(paths[0], 'Maps', 'Encoding', sub + '_' + roi + '_' + '+'.join(thisContrast) + '_P2Lalphas.nii.gz'))
    map2nifti(thisDS, calphas[1]).to_filename(
        os.path.join(paths[0], 'Maps', 'Encoding', sub + '_' + roi + '_' + '+'.join(thisContrast) + '_L2Palphas.nii.gz'))