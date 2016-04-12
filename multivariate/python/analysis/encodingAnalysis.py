#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
# this script represents just throwing pymvpa at the problem. doesn't work great, and I suspect it's
# because we're using an encoding model.
"""Encoding analysis: using only motion corrected and coregistered data. This analysis applies a Savitsky-Golay filter,
then...
todo: add ridge regression, maybe try crossmodal classification"""
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
    from mvpa2.datasets import eventrelated as er
    TR = np.median(np.diff(thisDS.sa.time_coords))
    idx = 0
    # events are loading wrong...
    theseEvents =beta_events[sub]
    events = []
    for i in np.arange(len(beta_events[sub])):
        for ev in theseEvents[i]:
            ev['chunks'] = thisDS.sa['chunks'].unique[i]
            ev['onset'] += TR*idx
            ev['targets'] = ev['condition']
            del ev['intensity']
            if ev['duration'] is not '0':
                events.append(ev)
        idx += np.sum(thisDS.sa['chunks'].value == thisDS.sa['chunks'].unique[i])

    evds = er.fit_event_hrf_model(thisDS, events, time_attr='time_coords',
                                  condition_attr=('targets', 'chunks'),
                                  design_kwargs={'add_regs': mc_params[sub], 'hrf_model': 'canonical'},
                                  return_model=True)

    fds = lmvpa.replacetargets(evds, contrasts, thisContrast)
    fds = fds[fds.sa.targets != '0'] # remove the probes

    # prototyping
    # from nipy.modalities.fmri.design_matrix import make_dmtx
    # make_dmtx(frametimes, paradigm=None, hrf_model='canonical', drift_model='cosine', hfcut=128, drift_order=1, fir_delays=[0], add_regs=None, add_reg_names=None, min_onset=-24)