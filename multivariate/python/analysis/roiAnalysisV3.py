#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
# this script represents just throwing pymvpa at the problem. doesn't work great, and I suspect it's
# because we're using an encoding model.
"""V3.0: using only motion corrected and coregistered data. This analysis applies a Savitsky-Golay filter,
runs a beta extraction for each trial, subs in the corresponding labels and classifies."""
# try doing whole brain classification with heavy feature selection
import sys
import numpy as np
# initialize stuff
if sys.platform == 'darwin':
    plat = 'mac'
    sys.path.append('/Users/njchiang/GitHub/LanguageMVPA/multivariate/python/utils')
    debug=True
else:
    plat = 'win'
    sys.path.append('D:\\GitHub\\LanguageMVPA\\multivariate\\python\\utils')
    debug=False
import lmvpautils as lmvpa
from string import Template
# plat = 'usb'
paths, subList, contrasts, maskList = lmvpa.initpaths(plat)
thisContrast = 'verb'
chance = .125
roi = 'grayMatter'
nVox = 1000  # number of voxels to select
nComp = 10  # number of SVD components
filterLen = 49  #for SG filter
filterOrd = 3  # for sg filter
if debug:
    subList = {'LMVPA005': subList['LMVPA005']}
# this would load all data at once. not memory efficient
# ds_all = lmvpa.loadsubdata(paths, subList, m=roi, c='trial_type')
# motion parameters for all subjects
mc_params = lmvpa.loadmotionparams(paths, subList)
# events for beta extraction
# add everything as a sample attribute
beta_events = lmvpa.loadevents(paths, subList)

from mvpa2.datasets.mri import map2nifti
from mvpa2.mappers.zscore import zscore
import SavGolFilter as sg
import mvpa2.datasets.eventrelated as er
from mvpa2.mappers.svd import SVDMapper
from mvpa2.datasets.base import StaticFeatureSelection, ChainMapper
from mvpa2.featsel.helpers import FixedNElementTailSelector
import mvpa2.featsel as fs
from mvpa2.clfs import svm
from mvpa2.clfs.warehouse import OneWayAnova, LDA
from mvpa2.clfs.meta import FeatureSelectionClassifier, MappedClassifier
from mvpa2.measures.base import CrossValidation
from mvpa2.generators.partition import NFoldPartitioner
from mvpa2.misc import errorfx
from mvpa2.clfs import gda
# import mvpa2.misc.plot as pt
#initialize classifier
# clf = lm.BayesianRidge()
# clf.fit(desX.matrix, thisDS.samples[:,1])
# clf.predict(desX.matrix)
# make class that takes skl regression models as input and wraps into mapper
cname = 'LinearSVM'
clf = svm.LinearCSVMC()
preproc = 'None'
# feature selection
# preproc = 'fsel-'+str(nVox)
# fsel=LDA()
fsel = OneWayAnova()
fselector = fs.helpers.FixedNElementTailSelector(nVox, tail='upper',
                                                 mode='select', sort=False)
# fselector = fs.helpers.FractionTailSelector(0.05, mode='select', tail='upper')
sbfs = fs.base.SensitivityBasedFeatureSelection(fsel, fselector,
                                                enable_ca=['sensitivities'])
fclf = FeatureSelectionClassifier(clf, sbfs)

# SVD
# preproc = 'SVD-' + str(nComp)
# svdmapper = SVDMapper()
# get_SVD_sliced = lambda x: ChainMapper([svdmapper, StaticFeatureSelection(x)])
# fclf = MappedClassifier(clf, get_SVD_sliced(slice(0, nComp)))

################################################
cv = CrossValidation(fclf,
                     NFoldPartitioner(attr='chunks'),
                     errorfx=errorfx.mean_match_accuracy)

lresults = []
presults = []
l2presults=[]
p2lresults=[]
labels = []
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

    #OLS beta extraction
    betas = 'LSA'
    evds = er.fit_event_hrf_model(rds, events, time_attr='time_coords',
                                  condition_attr=('trial_type',  'chunks'),
                                  design_kwargs={'add_regs': mc_params[sub], 'hrf_model': 'canonical'},
                                  return_model=True)

    fds = lmvpa.replacetargets(evds, contrasts, thisContrast)
    fds = fds[fds.targets != '0']
    fds = fds[fds.targets != 'rest']
    lidx = fds.chunks < fds.sa['chunks'].unique[len(fds.sa['chunks'].unique)/2]
    pidx = fds.chunks >= fds.sa['chunks'].unique[len(fds.sa['chunks'].unique) / 2]

    lfds = fds[lidx]
    pfds = fds[pidx]

    lres = cv(lfds)
    lresults.append(lres)
    print "language: " + str(np.mean(lres.samples))
    cv.untrain()
    pres = cv(pfds)
    presults.append(pres)
    print "pictures: " + str(np.mean(pres.samples))
    labels.append(sub)
    fclf.train(fds[lidx])
    l2presults.append(np.mean(fclf.predict(fds[pidx]) == fds[pidx].sa.targets))
    fclf.untrain()
    fclf.train(fds[pidx])
    p2lresults.append(np.mean(fclf.predict(fds[lidx]) == fds[lidx].sa.targets))


def plotsubs(lr, pr, l2p, p2l, c=None, title=None, bar_width=.2, opacity=.4, error_config={'ecolor': '0.3'}):
    import matplotlib.pyplot as plt
    # results is the concatenated output of cv across subjects... or something.
    f, (ax1, ax2) = plt.subplots(2, figsize=(12,6))
    index = np.arange(len(lr))
    lheights = []
    lerrbars = []
    pheights = []
    perrbars = []

    for i in lr:
        lheights.append(np.mean(i.samples))
        lerrbars.append(np.std(i.samples))
    lheights = np.array(lheights)
    lstd = np.array(lerrbars)

    rects1 = ax1.bar(index, lheights, bar_width,
                     alpha=opacity,
                     color='b',
                     yerr=lstd,
                     error_kw=error_config,
                     label='Language')
    for i in pr:
        pheights.append(np.mean(i.samples))
        perrbars.append(np.std(i.samples))
    pheights=np.array(pheights)
    pstd = np.array(perrbars)
    rects2 = ax1.bar(index+bar_width, pheights, bar_width,
                     alpha=opacity,
                     color='r',
                     yerr=pstd,
                     error_kw=error_config,
                     label='Picture')


    l2pheights=np.array(l2p)
    rects3 = ax1.bar(index+bar_width+bar_width, l2pheights, bar_width,
                     alpha=opacity,
                     color='g',
                     error_kw=error_config,
                     label='L2P')

    p2lheights=np.array(p2l)
    rects4 = ax1.bar(index+bar_width+bar_width+bar_width, p2lheights, bar_width,
                     alpha=opacity,
                     color='y',
                     error_kw=error_config,
                     label='P2L')

    if not c is None:
        ax1.plot((0, len(lr)), (c,c), 'k-', alpha=opacity/2, label='chance')
        ax2.plot((0, len(lr)), (c,c), 'k-', alpha=opacity/2, label='chance')
    ax1.set_xlabel('Subjects')
    ax1.set_ylabel('Accuracy')
    ax2.set_ylabel('Accuracy')
    if not desc is None:
        ax1.set_title(title)
    else:
        ax1.set_title('Classification Accuracies')

    ax1.legend()
    # ax2.legend()
    ax1.set_ylim(0, 1)
    ax2.set_ylim(0, 1)
    ax2.boxplot([lheights, pheights, l2pheights, p2lheights], labels=['Language', 'Pictures', 'L2P', 'P2L'])
    plt.show()


desc = Template('$con: SG filter($fl, $fo) Classifier($c) Preproc($p) Mask($m) Betas($b)')
title = desc.substitute(con=thisContrast, fl=str(filterLen), fo=str(filterOrd), c=cname, p=preproc, m=roi, b=betas)
plotsubs(lresults, presults, l2presults, p2lresults, c=chance, title=title)
