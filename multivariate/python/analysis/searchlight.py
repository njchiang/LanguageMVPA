#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
# this script represents just throwing pymvpa at the problem. doesn't work great, and I suspect it's
# because we're using an encoding model.
"""Encoding analysis: using only motion corrected and coregistered data. This analysis applies a Savitsky-Golay filter,
then...
todo: add ridge regression, maybe try crossmodal classification"""
# initialize stuff
import sys
import getopt

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
from sklearn import linear_model as lm
from mvpa2.clfs.skl import SKLLearnerAdapter
from mvpa2.clfs import svm
# debug = True

paths, subList, contrasts, maskList = lmvpa.initpaths(plat)
if debug:
    subList = {'LMVPA005': subList['LMVPA005']}

import numpy as np
import os
from mvpa2.datasets.mri import map2nifti
from mvpa2.mappers.zscore import zscore
import SavGolFilter as sg
import mvpa2.datasets.eventrelated as er
from mvpa2.measures.base import CrossValidation
from mvpa2.generators.partition import NFoldPartitioner
import searchlightutils as sl

# for sub in subList.keys():
def runsub(sub, thisContrast, r, dstype='raw', roi='grayMatter', filterLen=49, filterOrd=3, write=False):

    if dstype == 'raw':
        outdir='PyMVPA'
        print "working with raw data"
        thisSub = {sub: subList[sub]}
        dsdict = lmvpa.loadsubdata(paths, thisSub, m=roi)
        thisDS = dsdict[sub]
        mc_params = lmvpa.loadmotionparams(paths, thisSub)
        beta_events = lmvpa.loadevents(paths, thisSub)
        # savitsky golay filtering
        sg.sg_filter(thisDS, filterLen, filterOrd)
        # gallant group zscores before regression.

        # zscore w.r.t. rest trials
        # zscore(thisDS, param_est=('targets', ['rest']), chunks_attr='chunks')
        # zscore entire set. if done chunk-wise, there is no double-dipping (since we leave a chunk out at a time).
        zscore(thisDS, chunks_attr='chunks')
        print "beta extraction"
        ## BETA EXTRACTION ##
        rds, events = lmvpa.amendtimings(thisDS.copy(), beta_events[sub])
        evds = er.fit_event_hrf_model(rds, events, time_attr='time_coords',
                                      condition_attr=('trial_type', 'chunks'),
                                      design_kwargs={'add_regs': mc_params[sub], 'hrf_model': 'canonical'},
                                      return_model=True)

        fds = lmvpa.replacetargets(evds, contrasts, thisContrast)
        fds = fds[fds.targets != '0']
    else:
        outdir=os.path.join('LSS', dstype)
        print "loading betas"
        fds = lmvpa.loadsubbetas(paths, sub, btype=dstype, m=roi)
        fds.sa['targets'] = fds.sa[thisContrast]
        zscore(fds, chunks_attr='chunks')

    fds = lmvpa.sortds(fds)
    print "searchlights"
    ## initialize classifier
    clf = svm.LinearNuSVMC()
    cv = CrossValidation(clf, NFoldPartitioner())
    from mvpa2.measures.searchlight import sphere_searchlight
    cvSL = sphere_searchlight(cv, radius=r)


    # now I have betas per chunk. could just correlate the betas, or correlate the predictions for corresponding runs
    lidx = fds.chunks < fds.sa['chunks'].unique[len(fds.sa['chunks'].unique)/2]
    pidx = fds.chunks >= fds.sa['chunks'].unique[len(fds.sa['chunks'].unique) / 2]

    lres = sl.run_cv_sl(cvSL, fds[lidx].copy(deep=False))
    pres = sl.run_cv_sl(cvSL, fds[pidx].copy(deep=False))

    if write:
        from mvpa2.base import dataset
        map2nifti(fds, dataset.vstack([lres, pres])).\
            to_filename(os.path.join(
                        paths[0], 'Maps', outdir,
                        sub + '_' + roi + '_' + thisContrast + '_cvsl.nii.gz'))

    del lres, pres, cvSL

    cvSL = sphere_searchlight(cv, radius=r)
    crossSet = fds.copy()
    crossSet.chunks[lidx] = 1
    crossSet.chunks[pidx] = 2
    cres = sl.run_cv_sl(cvSL, crossSet.copy(deep=False))
    if write:
        map2nifti(fds, cres[0]).to_filename(
            os.path.join(paths[0], 'Maps', outdir,
                         sub + '_' + roi + '_' + (thisContrast) + '_P2L.nii.gz'))
        map2nifti(fds, cres[1]).to_filename(
            os.path.join(paths[0], 'Maps', outdir,
                         sub + '_' + roi + '_' + (thisContrast) + '_L2P.nii.gz'))



def main(argv):
    thisContrast = None
    debug = False
    write = False
    roi = 'grayMatter'
    dstype='raw'
    r = 4  # searchlight radius
    try:
        opts, args = getopt.getopt(argv, "dwhm:c:b:r", ["mfile=", "contrast=", "debug="])
    except getopt.GetoptError:
        print 'searchlight.py -m <maskfile> -c <contrast> '
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'searchlight.py -m <maskfile> -c <contrast>  '
            sys.exit()
        elif opt in ("-m", "--mask"):
            roi = arg
        # elif opt in ("-n", "--niter"):
        #     n = arg
        elif opt in ("-c", "--contrast"):
            thisContrast = arg
        elif opt in ("-d", "--debug"):
            print "debug mode"
            debug = True
        elif opt in ("-w", "--write"):
            write = True
        elif opt in ("-r", "--radius"):
            r = arg
        elif opt in ("-b", "--dstype"):
            dstype = arg

    if not thisContrast:
        print "not a valid contrast... exiting"
        sys.exit(1)

    paths, subList, contrasts, maskList = lmvpa.initpaths(plat)
    print "Mask: " + str(roi)
    print "Full Model: " + str(thisContrast)
    print "Searchlight Radius: " + str(r)
    print "Write results: " + str(write)
    sg_params = [49, 2]
    if debug:
        subList = {'LMVPA005': subList['LMVPA005']}

    for s in subList.keys():
        runsub(sub=s, thisContrast=thisContrast, r=r, dstype=dstype, write=write,
               filterLen=sg_params[0], filterOrd=sg_params[1], roi=roi)

if __name__ == "__main__":
    main(sys.argv[1:])