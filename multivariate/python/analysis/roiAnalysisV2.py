#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
# this script represents just throwing pymvpa at the problem. doesn't work great, and I suspect it's
# because we're using an encoding model.
"""V2.0: using betas extracted with FSL. this differs from V3 in that the preprocessing is different.
 Spatial smoothing: 5mm, hpf 100s, LSA extraction using FSL"""
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

import os

print "Initializing..."
# initialize paths
import lmvpautils as lmvpa
paths, subList, contrasts, maskList = lmvpa.initpaths(plat)
# nVox = 100
# initialize subjects, masks, contrast
roi = "left_IFG_operc"
con = "syntax"
chance = .25
dsType = "Lang"
nVox = 500
# dsType = "Pic"

cname = 'LinearSVM'
from mvpa2.clfs import svm
clf = svm.LinearCSVMC()
preproc = 'None'
# feature selection
# preproc = 'fsel-'+str(nVox)
# fsel=LDA()
from mvpa2.clfs.warehouse import OneWayAnova, LDA
fsel = OneWayAnova()
import mvpa2.featsel as fs
fselector = fs.helpers.FixedNElementTailSelector(nVox, tail='upper',
                                                 mode='select', sort=False)
# fselector = fs.helpers.FractionTailSelector(0.05, mode='select', tail='upper')
sbfs = fs.base.SensitivityBasedFeatureSelection(fsel, fselector,
                                                enable_ca=['sensitivities'])
from mvpa2.clfs.meta import FeatureSelectionClassifier, MappedClassifier
fclf = FeatureSelectionClassifier(clf, sbfs)

from mvpa2.measures.base import CrossValidation
from mvpa2.misc import errorfx
from mvpa2.generators.partition import NFoldPartitioner

cv = CrossValidation(fclf,
                     NFoldPartitioner(attr='chunks'),
                     errorfx=errorfx.mean_match_accuracy)

import numpy as np
from mvpa2.misc.io.base import SampleAttributes
cv_attr = SampleAttributes(os.path.join(paths[3], (con + "_attribute_labels.txt")))

from mvpa2.measures import rsa
dsm = rsa.PDist(square=True)
# searchlight
# import searchlightutils as sl
# from mvpa2.measures.searchlight import sphere_searchlight
# cvSL = sphere_searchlight(cv, radius=r)
# lres = sl.run_cv_sl(cvSL, fds[lidx].copy(deep=False))

lresults = []
presults = []
l2presults = []
p2lresults = []
rsaresults = []
for sub in subList.keys():
    betas = lmvpa.loadsubbetas(paths, sub, m=roi, a=cv_attr)
    rsaresults.append(dsm(betas))

    # should i zscore?
    lidx = np.arange(32)
    pidx = np.arange(32, 64)

    lres = cv(betas[lidx].copy())
    lresults.append(lres)
    print "language: " + str(np.mean(lres.samples))
    cv.untrain()
    pres = cv(betas[pidx].copy())
    presults.append(pres)
    print "pictures: " + str(np.mean(pres.samples))
    cv.untrain()

    fclf.train(betas[lidx].copy())
    l2presults.append(np.mean(fclf.predict(betas[pidx]) == betas[pidx].sa.targets))
    fclf.untrain()
    fclf.train(betas[pidx])
    p2lresults.append(np.mean(fclf.predict(betas[lidx]) == betas[lidx].sa.targets))


import matplotlib.pyplot as plt
def plotsubs(lr, pr, l2p, p2l, c=None, title=None, bar_width=.2, opacity=.4, error_config={'ecolor': '0.3'}):
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

from string import Template
desc = Template('$con: Classifier($c) Mask($m)')
title = desc.substitute(con=con,  c=cname, m=roi)
plotsubs(lresults, presults, l2presults, p2lresults, c=chance, title=title)

import rsautils as ru
ru.plot_mtx(ru.rankTransform(np.mean(np.dstack(rsaresults), axis=2)), betas.sa.targets, 'ROI pattern correlation distances')
