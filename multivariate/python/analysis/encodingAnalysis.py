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
import numpy as np
# plat = 'usb'
# debug = True
paths, subList, contrasts, maskList = lmvpa.initpaths(plat)
thisContrast = 'syntax'
roi = 'left_IFG_operc'
if debug:
    subList = {'LMVPA002': subList['LMVPA002']}

# load things in as trial type for easy regression, then swap out labels accordingly
# do we actually want to load all data simultaneously? or one brain at a time...
# ds_all = lmvpa.loadsubdata(paths, subList, m=roi, c='trial_type')
# motion parameters for all subjects
mc_params = lmvpa.loadmotionparams(paths, subList)
# events for beta extraction
beta_events = lmvpa.loadevents(paths, subList, c='trial_type')
beta_events = lmvpa.loadevents(paths, subList, c='verb')


######### for testing
# sub = 'LMVPA002'
# ds = ds_all[sub]
# from mvpa2.datasets.miscfx import summary
# print summary(ds)
# lmvpa.testsg(ds, 79, 3, 10, c='chunks')
###############
def runCVBootstrap(ds, X, part='chunks', nchunks=2, nboots=100, alphas=None):
    # runs cross validation on the chunks of the dataset (leave-one-out)
    import BootstrapRidgeMapper as bsr
    import numpy as np
    if alphas == None:
        alphas = np.logspace(0, 3, 20)
    # push design into source dataset
    # mds=[]
    res=[]
    from mvpa2.datasets import Dataset
    for c in np.arange(len(ds.sa[part].unique)):
        # need to combine regressors into one...
        trainidx =ds.sa[part].value != ds.sa[part].unique[c]
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


    from mvpa2.base import dataset
    return Dataset(np.vstack(res), sa={part: ds.sa[part].unique})


    # later this will loop
for sub in subList.keys():
    thisSub={sub: subList[sub]}
    ds = lmvpa.loadsubdata(paths, thisSub, m=roi, c='trial_type')
    thisDS=ds[sub]
    # savitsky golay filtering
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
    # desX = lmvpa.make_designmat(rds, events, time_attr='time_coords', condition_attr=['targets', 'chunks'],
    #                             design_kwargs={'add_regs': mc_params[sub], 'hrf_model': 'canonical'},
    #                             glmfit_kwargs=None, regr_attrs=None)



    # GLM
    # normal regression. doesn't use desX from above.
    import mvpa2.datasets.eventrelated as er
    evds = er.fit_event_hrf_model(rds, events, time_attr='time_coords',
                                  condition_attr=('targets', 'chunks'),
                                  design_kwargs={'add_regs': mc_params[sub], 'hrf_model': 'canonical'},
                                  return_model=True)

    # Ridge
    desX = lmvpa.make_designmat(rds, events, time_attr='time_coords', condition_attr=['targets'],
                                design_kwargs={'hrf_model': 'canonical', 'drift_model': 'blank'},
                                glmfit_kwargs=None, regr_attrs=None)

    import copy
    lidx=thisDS.chunks<thisDS.sa['chunks'].unique[len(thisDS.sa['chunks'].unique)/2]
    lds=copy.copy(desX)
    lds.matrix=lds.matrix[lidx]
    lres = runCVBootstrap(rds.copy()[lidx], lds)
    print 'language ' + str(np.mean(lres))
    pidx = thisDS.chunks >= thisDS.sa['chunks'].unique[len(thisDS.sa['chunks'].unique) / 2]
    pds = copy.copy(desX)
    pds.matrix = pds.matrix[pidx]
    pres = runCVBootstrap(rds.copy()[pidx], pds)
    print 'pictures: ' + str(np.mean(pres))

    # now can make this happen for cross-modal stuff...

# maybe add this stuff later.
# # some regressors might be corresponding not to original condition_attr
# # so let's separate them out
# regressor_names = model_params.sa[glm_condition_attr].value
# condition_regressors = np.array([v in glm_condition_attr_map.values()[0]
#                                  for v in regressor_names])
# assert (condition_regressors.dtype == np.bool)
# if not np.all(condition_regressors):
#     # some regressors do not correspond to conditions and would need
#     # to be taken into a separate dataset
#     model_params.a['add_regs'] = model_params[~condition_regressors]
#     # then we process the rest
#     model_params = model_params[condition_regressors]
#     regressor_names = model_params.sa[glm_condition_attr].value
#
# # now define proper condition sa's
# for con, con_map in glm_condition_attr_map.iteritems():
#     model_params.sa[con] = [con_map[v] for v in regressor_names]
# model_params.sa.pop(glm_condition_attr)  # remove generated one
# return model_params


