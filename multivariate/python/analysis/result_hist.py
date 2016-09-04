import matplotlib.pyplot as plt
import numpy as np
import sys
import os
# initialize stuff
if sys.platform == 'darwin':
    plat = 'usb'
    # plat = 'mac'
    sys.path.append('/Users/njchiang/GitHub/python-fmri-utils/utils')
    debug = True
else:
    plat = 'win'
    sys.path.append('D:\\GitHub\\LanguageMVPA\\multivariate\\python\\analysis')
    sys.path.append('D:\\GitHub\\python-fmri-utils\\utils')
    debug = False

import lmvpautils as lmvpa
# paths, subList, contrasts, maskList, fileMap = bu.initpaths(plat)
from mvpa2.datasets.mri import fmri_dataset


def result_hist(scan, masksize):
    paths, subList, contrasts, masks = lmvpa.initpaths(plat)
    bfn = os.path.join(paths[0], 'Maps', 'Encoding', scan + '.nii.gz')
    f, axarr = plt.subplots(len(masks), sharex=True)
    axarr[0].set_title('Model: ' + scan + ' | Mask size: ' + str(masksize))
    for i, m in enumerate(masks):
        mask = os.path.join(paths[0], 'data', s, 'masks', s + '_' + m + '.nii.gz')
        d = fmri_dataset(bfn, mask=mask)
        axarr[i].hist(np.mean(d[0].samples,0), bins=30)
        axarr[i].axvline(x=np.median(d.samples), color='r')
        axarr[i].set_ylabel(m)

masksize='8mm'
result_hist('AA01', 'AA01_3mm_grayMatter_teststim+testside_cvAlpha_prior_ridge', masksize)
result_hist('AA01', 'AA01_3mm_grayMatter_teststim+testside_cvAlpha_noprior_ridge', masksize)
result_hist('AA01', 'AA01_3mm_grayMatter_trial_type_cvAlpha_noprior_ridge', masksize)
result_hist('AA01', 'AA01_3mm_grayMatter_trial_type_cvAlpha_prior_ridge', masksize)


covarmat = (np.matrix('1, 0, .5, .5, 0, 0; \
                                    0, 1, 0, 0, 0, 0; \
                                    .5, 0, 1, -1, 0, 0; \
                                    .5, 0, -1, 1, 0, 0; \
                                    0, 0, 0, 0, 1, 0; \
                                    0, 0, 0, 0, 0, 1'))

fig, ax1 = plt.subplots()
mapp = ax1.imshow(covarmat, interpolation='nearest')

xtickNames = plt.setp(ax1, xticklabels=['', 'morph', 'fix', 'run', 'walk', 'left', 'right'])
plt.setp(xtickNames, rotation=45, fontsize=12)
ytickNames=plt.setp(ax1, yticklabels=['', 'morph', 'fix', 'run', 'walk', 'left', 'right'])
plt.setp(xtickNames, fontsize=12)
plt.colorbar(mapp)