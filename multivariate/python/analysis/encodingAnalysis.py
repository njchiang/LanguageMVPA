from os.path import join as pjoin
import sys
import numpy as np
sys.path.append('D:\\GitHub\\LanguageMVPA\\multivariate\\python\\utils')
sys.path.append('/Users/njchiang/GitHub/LanguageMVPA/multivariate/python/utils')
import lmvpautils as lmvpa
# initialize stuff
if sys.platform == 'darwin':
    plat = 'mac'
else:
    plat = 'win'

# plat = 'usb'
paths, subList, contrasts, maskList = lmvpa.initpaths(plat)

subjs = {'LMVPA003': subList['LMVPA003']}
# load things in as trial type for easy regression, then swap out labels accordingly
ds_all = lmvpa.loadsubdata(paths, subjs, m='left_IFG_operc', c='trial_type')
mc_params = lmvpa.loadmotionparams(paths, subList)
# events for beta extraction
beta_events = lmvpa.loadevents(paths, subjs, c='trial_type')

######### for testing
# sub = 'LMVPA002'
# ds = ds_all[sub]
# from mvpa2.datasets.miscfx import summary
# print summary(ds)
# lmvpa.testsg(ds, 79, 3, 10, c='chunks')
###############

    # later this will loop
for sub in subjs.keys():

    # savitsky golay filtering
    thisDS = ds_all[sub].copy()
    import SavGolFilter as sg
    sg.sg_filter(thisDS, 79, 3)


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

    TR = np.median(np.diff(ds.sa.time_coords))
    idx = 0
    events=[]
    for i in np.arange(len(beta_events[sub])):
        for ev in beta_events[sub][i]:
            ev['chunks'] = thisDS.sa['chunks'].unique[i]
            ev['onset'] += TR*idx
            ev['targets'] = ev['condition']
            del ev['intensity']
            events.append(ev)
        idx += np.sum(thisDS.sa['chunks'].value == thisDS.sa['chunks'].unique[i])

    # # regress, getting timing data from dataset
    # events2 = er.find_events(targets=thisDS.sa.targets, chunks=thisDS.sa.chunks)
    # # events = [ev for ev in events if ev['targets'] not in ['rest']]
    #
    # # convert onsets and durations into timestamps
    # i = 0
    # for ev in events2:
    #     ev['onset'] *= TR
    #     ev['duration'] *= TR


    # maybe add confound regressors. how to get trial by trial estimate?
    # or, use nipy to generate design matrix and feed into ridge regression
    # rds = thisDS.copy(deep=False)
    # rds = er.assign_conditionlabels(rds, events, noinfolabel='rest')

    evds = er.fit_event_hrf_model(thisDS, events, time_attr='time_coords',
                                  condition_attr=('targets', 'chunks'),
                                  design_kwargs={'add_regs': mc_params[sub], 'hrf_model': 'canonical'})

    # thisDS = lmvpa.replacetargets(thisDS, contrasts, 'syntax')
    fds=lmvpa.replacetargets(evds, contrasts, 'syntax')

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
