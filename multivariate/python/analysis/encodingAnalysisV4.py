#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""
Current working version of encoding analysis with L2 regularization (when no prior structure is given, this is simply ridge regression.
call from command line to flexibly change the mask and feature space. This way the file doesn't have to be continually modified.
"""
import os
from mvpa2.datasets.mri import map2nifti
from mvpa2.mappers.zscore import zscore
import copy as cp
import numpy as np
import logging
import getopt

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
import SavGolFilter as sg
import BootstrapRidge as bsr

paths, subList, contrasts, maskList = lmvpa.initpaths(plat)


def runsub(sub, thisContrast, thisContrastStr,
           filterLen, filterOrd,
           paramEst, chunklen, alphas=np.logspace(0, 3, 20), roi='grayMatter'):
    thisSub = {sub: subList[sub]}
    mc_params = lmvpa.loadmotionparams(paths, thisSub)
    beta_events = lmvpa.loadevents(paths, thisSub)
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
    # rds, events = lmvpa.amendtimings(thisDS.copy(), beta_events[sub])
    rds, events = lmvpa.amendtimings(thisDS.copy(), beta_events[sub], contrasts) # adding features

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
    nchunks = int(len(thisDS)*paramEst / chunklen)

    covarmat = None
    mus = None
    lwts, lalphas, lres = bsr.bootstrap_ridge(rds[lidx], ldes, chunklen=chunklen, nchunks=nchunks,
                                              cov0=covarmat, mu0=mus, part_attr='chunks', mode='test',
                                              alphas=alphas, single_alpha=True, normalpha=False,
                                              nboots=15, corrmin=.2, singcutoff=1e-10, joined=None,
                                              use_corr=True)
    pwts, palphas, pres = bsr.bootstrap_ridge(rds[pidx], pdes, chunklen=chunklen, nchunks=nchunks,
                                              part_attr='chunks', mode='test',
                                              alphas=alphas, single_alpha=True, normalpha=False,
                                              nboots=15, corrmin=.2, singcutoff=1e-10, joined=None,
                                              use_corr=True)
    print 'language ' + str(np.mean(lres))

    # pictures within
    print 'pictures: ' + str(np.mean(pres))

# need to change outstring
    from mvpa2.base import dataset
    map2nifti(thisDS, dataset.vstack([lres, pres])) \
        .to_filename(os.path.join(paths[0], 'Maps', 'Encoding', sub + '_' + roi + '_' + thisContrastStr +
                                  '_ridge_corrs.nii.gz'))
    map2nifti(thisDS, dataset.vstack([lwts, pwts])) \
        .to_filename(os.path.join(paths[0], 'Maps', 'Encoding', sub + '_' + roi + '_' + thisContrastStr +
                                  '_ridge_weights.nii.gz'))
    map2nifti(thisDS, dataset.vstack([lalphas, palphas])) \
        .to_filename(os.path.join(paths[0], 'Maps', 'Encoding', sub + '_' + roi + '_' + thisContrastStr +
                                  '_ridge_alphas.nii.gz'))

    del lres, pres, lwts, pwts, lalphas, palphas
    crossSet = thisDS.copy()
    crossSet.chunks[lidx] = 1
    crossSet.chunks[pidx] = 2
    cwts, calphas, cres = bsr.bootstrap_ridge(crossSet, des, chunklen=chunklen, nchunks=nchunks,
                                              part_attr='chunks', mode='test',
                                              alphas=alphas, single_alpha=True, normalpha=False,
                                              nboots=15, corrmin=.2, singcutoff=1e-10, joined=None,
                                              use_corr=True)
    print 'cross: ' + str(np.mean(cres))
    map2nifti(thisDS, cres[0]).to_filename(
        os.path.join(paths[0], 'Maps', 'Encoding', sub + '_' + roi + '_' + thisContrastStr + '_P2L_ridge_corr.nii.gz'))
    map2nifti(thisDS, cres[1]).to_filename(
        os.path.join(paths[0], 'Maps', 'Encoding', sub + '_' + roi + '_' + thisContrastStr + '_L2P_ridge_corr.nii.gz'))

    map2nifti(thisDS, cwts[0]).to_filename(
        os.path.join(paths[0], 'Maps', 'Encoding', sub + '_' + roi + '_' + thisContrastStr + '_P2L_ridge_weights.nii.gz'))
    map2nifti(thisDS, cwts[1]).to_filename(
        os.path.join(paths[0], 'Maps', 'Encoding', sub + '_' + roi + '_' + thisContrastStr + '_L2P_ridge_weights.nii.gz'))

    map2nifti(thisDS, calphas[0]).to_filename(
        os.path.join(paths[0], 'Maps', 'Encoding', sub + '_' + roi + '_' + thisContrastStr + '_P2L_ridge_alphas.nii.gz'))
    map2nifti(thisDS, calphas[1]).to_filename(
        os.path.join(paths[0], 'Maps', 'Encoding', sub + '_' + roi + '_' + thisContrastStr + '_L2P_ridge_alphas.nii.gz'))

    del cres, cwts, calphas


def main(argv):
    roi = 'grayMatter'
    thisContrast = []

    paths, subList, contrasts, maskList = lmvpa.initpaths(plat)

    try:
      opts, args = getopt.getopt(argv, "hm:c:", ["mfile=", "contrast="])
    except getopt.GetoptError:
      print 'encodingAnalysisV4.py -r <roi> -m <maskfile> -c <contrast> '
      sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'encodingAnalysisV4.py -m <mask> -c <contrast> '
            sys.exit()
        elif opt in ("-m", "--mask"):
            roi = arg
        elif opt in ("-c", "--contrast"):
            thisContrast = arg.split(',')

    debug=False
    if not thisContrast:
        print "not a valid contrast... exiting"
        sys.exit(1)

    thisContrastStr = '+'.join(thisContrast)
    print(thisContrastStr)
    if 'word2vec' in thisContrast:
        thisContrast.remove('word2vec')
        for i in np.arange(0, 300):
            thisContrast.append('word2vec' + str(i))

    sg_params = [49, 2]
    chunklen = 12  # this reflects the length of a complete trial
    paramEst = .2  # this much data to be held out for ridge regression parameter estimation
    if debug:
        subList = {'LMVPA005': subList['LMVPA005']}

    # ds_all = lmvpa.loadsubdata(paths, subList, m=roi, c='trial_type')
    # motion parameters for all subjects
    # events for beta extraction
    # add everything as a sample attribute

    logging.basicConfig(level=logging.DEBUG)
    alphas = np.logspace(-3, 3, 50)
    for s in subList.keys():
        runsub(sub=s, thisContrast=thisContrast, thisContrastStr=thisContrastStr,
               filterLen=sg_params[0], filterOrd=sg_params[1],
               paramEst=paramEst, chunklen=chunklen, alphas=alphas, roi='grayMatter')

if __name__ == "__main__":
    main(sys.argv[1:])
