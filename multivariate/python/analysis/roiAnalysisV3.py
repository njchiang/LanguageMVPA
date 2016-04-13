#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
# this script represents just throwing pymvpa at the problem. doesn't work great, and I suspect it's
# because we're using an encoding model.
"""V3.0: using only motion corrected and coregistered data. This analysis applies a Savitsky-Golay filter,
runs a beta extraction for each trial, subs in the corresponding labels and classifies."""
import sys
import numpy as np
# initialize stuff
if sys.platform == 'darwin':
    plat = 'mac'
    sys.path.append('/Users/njchiang/GitHub/LanguageMVPA/multivariate/python/utils')
else:
    plat = 'win'
    sys.path.append('D:\\GitHub\\LanguageMVPA\\multivariate\\python\\utils')

import lmvpautils as lmvpa

# plat = 'usb'
paths, subList, contrasts, maskList = lmvpa.initpaths(plat)
thisContrast='syntax'
roi = 'left_IFG_operc'
# subList = {'LMVPA001': subList['LMVPA001']}
# load things in as trial type for easy regression, then swap out labels accordingly
ds_all = lmvpa.loadsubdata(paths, subList, m=roi, c='trial_type')
# motion parameters for all subjects
mc_params = lmvpa.loadmotionparams(paths, subList)
# events for beta extraction
beta_events = lmvpa.loadevents(paths, subList, c='trial_type')

######### for testing
# sub = 'LMVPA001'
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

    desX = lmvpa.make_designmat(rds, events, time_attr='time_coords', condition_attr=['targets', 'chunks'],
                                design_kwargs={'add_regs': mc_params[sub], 'hrf_model': 'canonical'},
                                glmfit_kwargs=None, regr_attrs=None)

    import mvpa2.datasets.eventrelated as er

    evds = er.fit_event_hrf_model(rds, events, time_attr='time_coords',
                                  condition_attr=('targets', 'chunks'),
                                  design_kwargs={'add_regs': mc_params[sub], 'hrf_model': 'canonical'},
                                  return_model=True)

    # # # regress, getting timing data from dataset
    # events2 = er.find_events(targets=thisDS.sa.targets, chunks=thisDS.sa.chunks)
    # # # events = [ev for ev in events if ev['targets'] not in ['rest']]
    # #
    # # # convert onsets and durations into timestamps
    # for ev in events2:
    #     ev['onset'] *= TR
    #     ev['duration'] *= TR
    #
    # events2 = [ev for ev in events2 if ev['targets'] != 'rest']
    # evds2 = er.fit_event_hrf_model(thisDS, events2, time_attr='time_coords',
    #                               condition_attr=('targets', 'chunks'),
    #                               design_kwargs={'add_regs': mc_params[sub], 'hrf_model': 'canonical'},
    #                               return_model=True)
    # or, use nipy to generate design matrix and feed into ridge regression

    fds = lmvpa.replacetargets(evds, contrasts, thisContrast)

    fds = fds[fds.sa.targets != '0']

    # split data into language and pictures
    from mvpa2.base import dataset
    lfds = [fds[fds.sa['chunks'].value == fds.sa['chunks'].unique[i]]
            for i in np.arange(len(fds.sa['chunks'].unique)/2)]
    lfds = dataset.vstack(lfds)

    pfds = [fds[fds.sa['chunks'].value == fds.sa['chunks'].unique[i]]
            for i in np.arange(len(fds.sa['chunks'].unique) / 2, len(fds.sa['chunks'].unique))]
    pfds = dataset.vstack(pfds)


    from mvpa2.clfs import svm
    clf = svm.LinearCSVMC()
    from mvpa2.measures.base import CrossValidation
    from mvpa2.generators.partition import NFoldPartitioner
    from mvpa2.misc import errorfx
    cv = CrossValidation(clf,
                         NFoldPartitioner(attr='chunks'),
                         errorfx=errorfx.mean_match_accuracy)
    cv.untrain()
    lres = cv(lfds)
    print "language: " + str(np.mean(lres.samples))
    cv.untrain()
    pres = cv(pfds)
    print "pictures: " + str(np.mean(pres.samples))
