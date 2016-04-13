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
debug = True
paths, subList, contrasts, maskList = lmvpa.initpaths(plat)
thisContrast = 'syntax'
roi = 'left_IFG_operc'
if debug:
    subList = {'LMVPA002': subList['LMVPA002']}

# load things in as trial type for easy regression, then swap out labels accordingly
# do we actually want to load all data simultaneously? or one brain at a time...
ds_all = lmvpa.loadsubdata(paths, subList, m=roi, c='trial_type')
# motion parameters for all subjects
mc_params = lmvpa.loadmotionparams(paths, subList)
# events for beta extraction
beta_events = lmvpa.loadevents(paths, subList, c='trial_type')

######### for testing
sub = 'LMVPA002'
# ds = ds_all[sub]
# from mvpa2.datasets.miscfx import summary
# print summary(ds)
# lmvpa.testsg(ds, 79, 3, 10, c='chunks')
###############

    # later this will loop
for sub in subList.keys():
    # savitsky golay filtering
    thisDS = ds_all[sub].copy()
    import SavGolFilter as sg
    sg.sg_filter(thisDS, 49, 3)
    # gallant group zscores before regression.
    # update events data
    from mvpa2.mappers.zscore import zscore
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
    desX = lmvpa.make_designmat(rds, events, time_attr='time_coords', condition_attr=['targets', 'chunks'],
                                design_kwargs={'add_regs': mc_params[sub], 'hrf_model': 'canonical'},
                                glmfit_kwargs=None, regr_attrs=None)

    # import ridge
    # import numpy as np
    # # estimate alpha via bootstrap
    # """    wt : array_like, shape (N, M)
    #     Regression weights for N features and M responses.
    # corrs : array_like, shape (M,)
    #     Validation set correlations. Predicted responses for the validation set are obtained using the regression
    #     weights: pred = np.dot(Pstim, wt), and then the correlation between each predicted response and each
    #     column in Presp is found.
    # alphas : array_like, shape (M,)
    #     The regularization coefficient (alpha) selected for each voxel using bootstrap cross-validation.
    # bootstrap_corrs : array_like, shape (A, M, B)
    #     Correlation between predicted and actual responses on randomly held out portions of the training set,
    #     for each of A alphas, M voxels, and B bootstrap samples.
    # valinds : array_like, shape (TH, B)
    #     The indices of the training data that were used as "validation" for each bootstrap sample.
    # """
    # wt, corrs, valphas, allRcorrs, valinds = ridge.bootstrap_ridge(Rstim=desX.matrix[rds.sa['chunks'].value==1],
    #                                                                Rresp=rds.samples[rds.sa['chunks'].value==1],
    #                                                                Pstim=desX.matrix[rds.sa['chunks'].value==2],
    #                                                                Presp=rds.samples[rds.sa['chunks'].value==2],
    #                                                                alphas=np.logspace(0, 3, 20),
    #                                                                nboots=100,
    #                                                                chunklen=4,
    #                                                                nchunks=4,
    #                                                                corrmin=0.2,
    #                                                                joined=None,
    #                                                                singcutoff=1e-10,
    #                                                                normalpha=False,
    #                                                                single_alpha=True,
    #                                                                use_corr=True)

    import BootstrapRidge as bsr
    import numpy as np
    import copy
    lidx = rds.chunks <= 2
    ldesX = copy.copy(desX)
    ldesX.matrix = ldesX.matrix[lidx]

    lridgeRes = bsr.bootstrap_ridge(rds[lidx], des=ldesX, alphas=np.logspace(0,3,20), nboots=100, chunklen=4, nchunks=4, single_alpha=True)
    lridgeResSummary = np.mean(lridgeRes.samples, axis=0)
    pdesX = copy.copy(desX)
    pdesX.matrix = pdesX.matrix[pidx]
    # find the correlations using bootstrapped alpha



    # make designmat. regress for each run and try to predict the other run
    # need to regress now...

    # normal regression. doesn't use desX from above.
    # make ridge GLM mapper, modify statsmodels_glm _fit_model
    import mvpa2.datasets.eventrelated as er
    evds = er.fit_event_hrf_model(rds, events, time_attr='time_coords',
                                  condition_attr=('targets', 'chunks'),
                                  design_kwargs={'add_regs': mc_params[sub], 'hrf_model': 'canonical'},
                                  return_model=True)

    fds = lmvpa.replacetargets(evds, contrasts, thisContrast)
    fds = fds[fds.sa.targets != '0'] # remove the probes


