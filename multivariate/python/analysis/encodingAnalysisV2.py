SKLLearnerAdapter(DecisionTreeRegressor(max_depth=2))


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
plat = 'usb'
# debug = True
paths, subList, contrasts, maskList = lmvpa.initpaths(plat)
thisContrast = 'syntax'
roi = 'grayMatter'
# if debug:
#     subList = {'LMVPA002': subList['LMVPA002']}

# load things in as trial type for easy regression, then swap out labels accordingly
# doing it one at a time is more memory efficient.
# ds_all = lmvpa.loadsubdata(paths, subList, m=roi, c='trial_type')
# motion parameters for all subjects
mc_params = lmvpa.loadmotionparams(paths, subList)
# events for beta extraction
beta_events = lmvpa.loadevents(paths, subList, c=thisContrast)

import copy
import os
from mvpa2.datasets.mri import map2nifti
from mvpa2.mappers.zscore import zscore
import SavGolFilter as sg
import mvpa2.datasets.eventrelated as er
for sub in subList.keys():
    thisSub={sub: subList[sub]}
    dsdict = lmvpa.loadsubdata(paths, thisSub, m=roi, c=thisContrast)
    thisDS=dsdict[sub]
    # savitsky golay filtering

    sg.sg_filter(thisDS, 49, 3)
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
    # now: make timing files for each feature that encompass every timepoint but with different intensities
    # desX = lmvpa.make_designmat(rds, events, time_attr='time_coords', condition_attr=['targets', 'chunks'],
    #                             design_kwargs={'add_regs': mc_params[sub], 'hrf_model': 'canonical'},
    #                             glmfit_kwargs=None, regr_attrs=None)



    # GLM
    # normal regression. doesn't use desX from above.

    evds = er.fit_event_hrf_model(rds, events, time_attr='time_coords',
                                  condition_attr=('targets', 'chunks'),
                                  design_kwargs={'add_regs': mc_params[sub], 'hrf_model': 'canonical'},
                                  return_model=True)

    # Ridge
    desX = lmvpa.make_designmat(rds, events, time_attr='time_coords', condition_attr=['targets'],
                                design_kwargs={'hrf_model': 'canonical', 'drift_model': 'blank'},
                                glmfit_kwargs=None, regr_attrs=None)


    from sklearn.linear_model import BayesianRidge
    from mvpa2.clfs.skl.base import SKLLearnerAdapter
    br = BayesianRidge(alpha_1=1e-06, alpha_2=1e-06, compute_score=False, copy_X=True,
       fit_intercept=True, lambda_1=1e-06, lambda_2=1e-06, n_iter=300,
       normalize=False, tol=0.001, verbose=False)
    SKLLearnerAdapter(br) # turn this into a mapper
    clf = BayesianRidge(compute_score=True)
    clf.fit(X, y)

    lidx=thisDS.chunks<thisDS.sa['chunks'].unique[len(thisDS.sa['chunks'].unique)/2]
    lds=copy.copy(desX)
    lds.matrix=lds.matrix[lidx]
    lres = runCVBootstrap(rds.copy()[lidx], lds)
    print 'language ' + str(np.mean(lres))
    map2nifti(thisDS, np.mean(lres, axis=0)).to_filename(os.path.join(paths[0], 'Maps', 'Encoding', sub + '_' + thisContrast + '_Lridge.nii.gz'))
    pidx = thisDS.chunks >= thisDS.sa['chunks'].unique[len(thisDS.sa['chunks'].unique) / 2]
    pds = copy.copy(desX)
    pds.matrix = pds.matrix[pidx]
    pres = runCVBootstrap(rds.copy()[pidx], pds)
    print 'pictures: ' + str(np.mean(pres))
    map2nifti(thisDS, np.mean(pres, axis=0)).to_filename(os.path.join(paths[0], 'Maps', 'Encoding', sub + '_' + thisContrast + '_Pridge.nii.gz'))
