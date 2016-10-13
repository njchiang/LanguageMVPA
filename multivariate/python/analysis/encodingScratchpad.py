import os
from mvpa2.datasets.mri import map2nifti
from mvpa2.mappers.zscore import zscore
import copy as cp
import numpy as np
import logging
import getopt
from nipy.modalities.fmri.design_matrix import make_dmtx
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
import matplotlib.pyplot as plt

paths, subList, contrasts, maskList = lmvpa.initpaths(plat)

sub = 'LMVPA005'

thisContrast = ['ap', 'cr']
thisContrastStr = 'ap+cr'
roi = 'grayMatter'
filterLen = 49
filterOrd = 2
chunklen = 30
paramEst = .25
thisSub = {sub: subList[sub]}
mc_params = lmvpa.loadmotionparams(paths, thisSub)
beta_events = lmvpa.loadevents(paths, thisSub)

###
dsdict = lmvpa.loadsubdata(paths, thisSub, m=roi, c='trial_type')
thisDS = dsdict[sub]
iVox = np.random.randint(0, thisDS.shape[1])
plt.plot(thisDS[:, iVox]-np.mean(thisDS[:, iVox]), label='raw')
sg.sg_filter(thisDS, filterLen, filterOrd)
plt.plot(thisDS[:, iVox], label='S-G')
zscore(thisDS, chunks_attr='chunks')
plt.plot(thisDS[:, iVox], label='Z')
plt.title('Voxel ' + str(iVox))
plt.legend()
###

rds, events = lmvpa.amendtimings(thisDS.copy(), beta_events[sub], contrasts)  # adding features

# we can model out motion and just not use those betas.
# Ridge
if isinstance(thisContrast, basestring):
    thisContrast = [thisContrast]
# instead of binarizing each one, make them parametric
desX, rds = lmvpa.make_designmat(rds, events, time_attr='time_coords', condition_attr=thisContrast,
                                 design_kwargs={'hrf_model': 'canonical', 'drift_model': 'blank'},
                                 regr_attrs=None)
desX, rds = lmvpa.make_designmat(rds, events, time_attr='time_coords', condition_attr=thisContrast,
                                 design_kwargs={'hrf_model': 'fir', 'fir_delays': [0, 1, 2, 3], 'drift_model': 'blank'},
                                 regr_attrs=None)
# want to collapse ap and cr, but have anim separate
desX['motion'] = make_dmtx(rds.sa['time_coords'].value, paradigm=None, add_regs=mc_params[sub], drift_model='blank')

des = lmvpa.make_parammat(cp.copy(desX), zscore=True)

