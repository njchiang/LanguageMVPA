#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""
Current working version of encoding analysis with L2 regularization
this takes a single alpha and compares each subject against a series of random design matrices
run ridge several times on random matrices and take average and save variance, and subtract actual corr from it
(when no prior structure is given, this is simply ridge regression.
call from command line to flexibly change the mask and feature space. This way the file doesn't have to be continually modified.
# so we'll use a good alpha and add mc params and test correlations for subsets of the model. use ridge to get everything
probably should normalize design matrix
# ignoring random condition in calculation for these betas. taken as the average of the following:
word2vec
probe+word2vec
ap+cr+word2vec
probe+ap+cr+word2vec
Language: 0.79061, 63.57399,
Pictures:  2.693528,  34.62643
best alpha overall:
best alpha language: 1.676833, 60.59053
best alpha pictures: 3.55648, 3.511738
"""
import os
from mvpa2.datasets.mri import map2nifti
from mvpa2.mappers.zscore import zscore
import copy as cp
import numpy as np
import logging
import getopt
from mvpa2.base import dataset
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


def runsub(sub, thisContrast, thisContrastStr, testContrast,
           filterLen, filterOrd, write=False, debug=False,
           alphas=1, roi='grayMatter'):
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

    # get timing data from timing files
    rds, events = lmvpa.amendtimings(thisDS.copy(), beta_events[sub], contrasts)  # adding features

    # we can model out motion and just not use those betas.
    # Ridge
    if isinstance(thisContrast, basestring):
        thisContrast = [thisContrast]
    desX, rds = lmvpa.make_designmat(rds, events, time_attr='time_coords', condition_attr=thisContrast,
                                     design_kwargs={'hrf_model': 'canonical', 'drift_model': 'blank'},
                                     regr_attrs=None)
    # 'add_regs': mc_params[sub]

    desX['motion'] = make_dmtx(rds.sa['time_coords'].value, paradigm=None, add_regs=mc_params[sub], drift_model='blank')

    des = lmvpa.make_parammat(desX, hrf='canonical', zscore=True)

    # split by language and pictures
    lidx = thisDS.chunks < thisDS.sa['chunks'].unique[len(thisDS.sa['chunks'].unique) / 2]
    pidx = thisDS.chunks >= thisDS.sa['chunks'].unique[len(thisDS.sa['chunks'].unique) / 2]
    ldes = cp.copy(des)
    pdes = cp.copy(des)

    ldes.matrix = ldes.matrix[lidx]
    pdes.matrix = pdes.matrix[pidx]

    covarmat = None
    mus = None
    lwts, _, lres, lceil = bsr.bootstrap_ridge(ds=rds[lidx], des=ldes, chunklen=1, nchunks=1,
                                               cov0=None, mu0=None, part_attr='chunks', mode='test',
                                               alphas=[alphas[0]], single_alpha=True, normalpha=False,
                                               nboots=1, corrmin=.2, singcutoff=1e-10, joined=None,
                                               use_corr=True)
    print 'language ' + str(np.mean(lres))

    pwts, _, pres, pceil = bsr.bootstrap_ridge(ds=rds[pidx], des=pdes, chunklen=1, nchunks=1,
                                               cov0=None, mu0=None, part_attr='chunks', mode='test',
                                               alphas=[alphas[1]], single_alpha=True, normalpha=False,
                                               nboots=1, corrmin=.2, singcutoff=1e-10, joined=None,
                                               use_corr=True)

    # pictures within
    print 'pictures: ' + str(np.mean(pres))
    if write:
        map2nifti(thisDS, dataset.vstack([lres, pres])) \
            .to_filename(os.path.join(paths[0], 'Maps', 'Encoding', sub + '_' + roi + '_' + thisContrastStr +
                                      '_ridge_la_' + str(alphas[0]) + '_pa_' + str(alphas[1]) + '_corrs.nii.gz'))
        map2nifti(thisDS, dataset.vstack([lwts, pwts])) \
            .to_filename(os.path.join(paths[0], 'Maps', 'Encoding', sub + '_' + roi + '_' + thisContrastStr +
                                      '_ridge_la_' + str(alphas[0]) + '_pa_' + str(alphas[1]) + '_wts.nii.gz'))
        map2nifti(thisDS, dataset.vstack([lceil, pceil])) \
            .to_filename(os.path.join(paths[0], 'Maps', 'Encoding', sub + '_' + roi + '_' + thisContrastStr +
                                      '_ridge_la_' + str(alphas[0]) + '_pa_' + str(alphas[1]) + '_ceiling.nii.gz'))

    for t in testContrast:
        tstr = '+'.join(t)
        lcorr = lmvpa.testmodel(wts=lwts, des=ldes, ds=rds[lidx], tc=cp.copy(t), use_corr=True)
        pcorr = lmvpa.testmodel(wts=pwts, des=pdes, ds=rds[pidx], tc=cp.copy(t), use_corr=True)
        if write:
            map2nifti(thisDS, dataset.vstack([lcorr, pcorr])) \
            .to_filename(os.path.join(paths[0], 'Maps', 'Encoding', sub + '_' + roi + '_' + tstr +
                                      '_ridge_la_' + str(alphas[0]) + '_pa_' + str(alphas[1]) + '_test_corrs.nii.gz'))

    del lres, pres, lwts, pwts, lceil, pceil
    crossSet = thisDS.copy()
    crossSet.chunks[lidx] = 1
    crossSet.chunks[pidx] = 2
    # cwts, cres, cceil = bsr.ridge(rds[pidx], pdes, mu0=mus, cov0=covarmat,
    #                                             part_attr='chunks', mode='test', alphas=alphas[0], single_alpha=True,
    #                                             normalpha=False, corrmin=.2, singcutoff=1e-10, joined=None,
    #                                             use_corr=True)
    cwts, _, cres, cceil = bsr.bootstrap_ridge(ds=crossSet, des=des, chunklen=1, nchunks=1,
                                               cov0=None, mu0=None, part_attr='chunks', mode='test',
                                               alphas=[alphas[2]], single_alpha=True, normalpha=False,
                                               nboots=1, corrmin=.2, singcutoff=1e-10, joined=None,
                                               use_corr=True)
    for t in testContrast:
        tstr = '+'.join(t)
        ccorr = lmvpa.testmodel(wts=cwts, des=des, ds=crossSet, tc=cp.copy(t), use_corr=True)
        if write:
            map2nifti(thisDS, ccorr[0]) \
            .to_filename(os.path.join(paths[0], 'Maps', 'Encoding', sub + '_' + roi + '_' + tstr +
                                      '_P2L_ridge_alpha_' + str(alphas[2]) + '_test_corr.nii.gz'))
            map2nifti(thisDS, ccorr[1]) \
            .to_filename(os.path.join(paths[0], 'Maps', 'Encoding', sub + '_' + roi + '_' + tstr +
                                      '_L2P_ridge_alpha_' + str(alphas[2]) + '_test_corr.nii.gz'))
    print 'cross: ' + str(np.mean(cres))
    if write:
        map2nifti(thisDS, cres[0]).to_filename(
        os.path.join(paths[0], 'Maps', 'Encoding', sub + '_' + roi + '_' + thisContrastStr +
                     '_P2L_ridge_alpha_' + str(alphas[2]) + '_corr.nii.gz'))
        map2nifti(thisDS, cres[1]).to_filename(
            os.path.join(paths[0], 'Maps', 'Encoding', sub + '_' + roi + '_' + thisContrastStr +
                         '_L2P_ridge_alpha_' + str(alphas[2]) + '_corr.nii.gz'))

        map2nifti(thisDS, cwts[cwts.chunks==1]).to_filename(
            os.path.join(paths[0], 'Maps', 'Encoding', sub + '_' + roi + '_' + thisContrastStr +
                         '_P2L_ridge_alpha_' + str(alphas[2]) + '_wts.nii.gz'))
        map2nifti(thisDS, cwts[cwts.chunks==2]).to_filename(
            os.path.join(paths[0], 'Maps', 'Encoding', sub + '_' + roi + '_' + thisContrastStr +
                         '_L2P_ridge_alpha' + str(alphas[2]) + '_wts.nii.gz'))

        map2nifti(thisDS, cceil[0]).to_filename(
            os.path.join(paths[0], 'Maps', 'Encoding', sub + '_' + roi + '_' + thisContrastStr +
                         '_P2L_ridge_alpha_' + str(alphas[2]) + '_ceiling.nii.gz'))
        map2nifti(thisDS, cceil[1]).to_filename(
            os.path.join(paths[0], 'Maps', 'Encoding', sub + '_' + roi + '_' + thisContrastStr +
                         '_L2P_ridge_alpha_' + str(alphas[2]) + '_ceiling.nii.gz'))
    del cres, cwts, cceil


def main(argv):
    roi = 'grayMatter'
    thisContrast = []
    testContrast = []
    debug=False
    write=False
    n = 100
    alphas = [1, 2.69, 2.62] # 1.68, 3.56]

    try:
        opts, args = getopt.getopt(argv, "dwhm:c:t:a", ["mfile=", "contrast=", "test=", "alpha=", "debug="])
    except getopt.GetoptError:
        print 'encodingAnalysisV5.py -m <maskfile> -c <contrast> -t <testcontrasts> -a <alphas>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'encodingAnalysisV5.py -m <maskfile> -c <contrast> -t <testcontrasts> -a <alphas> '
            sys.exit()
        elif opt in ("-m", "--mask"):
            roi = arg
        # elif opt in ("-n", "--niter"):
        #     n = arg
        elif opt in ("-c", "--contrast"):
            thisContrast = arg.split(',')
        elif opt in ("-t", "--test"):
            tmp = arg.split(',')
            testContrast.append(tmp)
        elif opt in ("-a", "--alpha"):
            tmp = arg.split(',')
            alphas = [float(i) for i in tmp]
        elif opt in ("-d", "--debug"):
            print "debug mode"
            debug = True
        elif opt in ("-w", "--write"):
            write = True

    if not thisContrast:
        print "not a valid contrast... exiting"
        sys.exit(1)

    paths, subList, contrasts, maskList = lmvpa.initpaths(plat)
    print "Mask: " + str(roi)
    print "Full Model: " + str(thisContrast)
    print "Model test: " + str(testContrast)
    print "Alphas: " + str(alphas)

    thisContrastStr = '+'.join(thisContrast)
    print(thisContrastStr)
    if 'word2vec' in thisContrast:
        thisContrast.remove('word2vec')
        for i in np.arange(0, 300):
            thisContrast.append('word2vec' + str(i))

    if 'random' in thisContrast:
        thisContrast.remove('random')
        for i in np.arange(0, 302):
            thisContrast.append('random' + str(i))

    sg_params = [49, 2]
    if debug:
        subList = {'LMVPA005': subList['LMVPA005']}

    logging.basicConfig(level=logging.DEBUG)
    for s in subList.keys():
        runsub(sub=s, thisContrast=thisContrast, thisContrastStr=thisContrastStr, debug=debug, write=write,
               testContrast=testContrast, filterLen=sg_params[0], filterOrd=sg_params[1],
               alphas=alphas, roi=roi)

if __name__ == "__main__":
    main(sys.argv[1:])
