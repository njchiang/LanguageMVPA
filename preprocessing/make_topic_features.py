
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
debug = True
thisContrast = ['anim', 'verb', 'ap', 'cr']
# thisContrast = ['verb', 'ap', 'cr']
# thisContrast = ['ap', 'cr']
# thisContrast = ['verb']
# thisContrast = ['anim', 'verb']
# thisContrast = ['anim', 'ap', 'cr']
# thisContrast = ['anim']

roi = 'grayMatter'
filterLen = 49
filterOrd = 3
chunklen=10 # this reflects the length of a complete trial
paramEst=.2 #this much data to be held out for ridge regression parameter estimation
paths, subList, contrasts, maskList = lmvpa.initpaths(plat)
if debug:
    subList = {'LMVPA005': subList['LMVPA005']}

# ds_all = lmvpa.loadsubdata(paths, subList, m=roi, c='trial_type')
# motion parameters for all subjects
mc_params = lmvpa.loadmotionparams(paths, subList)
# events for beta extraction
# add everything as a sample attribute
beta_events = lmvpa.loadevents(paths, subList)

# import BootstrapRidgeMapper as bsr
import numpy as np
import os
import matplotlib.pyplot as plt
import sklearn.decomposition as dcmp

alphas = np.logspace(0, 3, 20)

words = ['gravy', 'trainer', 'candle', 'bus', 'tree', 'goalie', 'singer', 'woman',
         'touched', 'stretched', 'lit', 'hit', 'crushed', 'kicked', 'kissed',
         'consoled', 'broccoli', 'dancer', 'match', 'truck', 'car', 'referee',
         'guitarist', 'man']

featureMat = np.loadtxt(os.path.join(paths[3], 'topics300_wordVecs.txt'), skiprows=1)[:, 1:].T
# normalize each row
featureMat = featureMat / np.sum(featureMat,1)[:, np.newaxis]
verbFeatureMat = featureMat[8:16,:]
verbs = words[8:16]
idx = np.array([7, 4, 3, 5, 6, 2, 1, 0])
verbsSort = [verbs[i] for i in idx]


f, axarr = plt.subplots(2, 3)

# manual feature reduction
thisFeatureMat = verbFeatureMat[:, np.sum(verbFeatureMat,0) > 0]
axarr[0,0].imshow(featureMat, aspect='auto', interpolation='nearest')

binFeatureMat = thisFeatureMat > 0
uniqueFeatureMat = thisFeatureMat[:, np.sum(binFeatureMat,0)>1]
axarr[0,1].imshow(uniqueFeatureMat, aspect='auto', interpolation='nearest')
im = axarr[0,2].imshow(np.cov(uniqueFeatureMat.T), aspect='equal', interpolation='nearest')

#PCs
thr=.95
p = dcmp.PCA()
p.fit(thisFeatureMat)
ncomponents = np.where(np.cumsum(p.explained_variance_ratio_) > thr)[0][0]+1 #because 0 indexing
del p
p = dcmp.PCA(n_components=ncomponents)
pcaFeatureMat = p.fit_transform(thisFeatureMat)
axarr[1,1].imshow(pcaFeatureMat, aspect='auto', interpolation='nearest')
im2 = axarr[1,2].imshow(np.cov(pcaFeatureMat.T), aspect='equal', interpolation='nearest')

umat = []
pmat = []
for i in np.arange(uniqueFeatureMat.shape[0]):
    umat.append(np.tile(uniqueFeatureMat[i,:], (4,1)))
    pmat.append(np.tile(pcaFeatureMat[i,:], (4,1)))
uniqueMat = np.tile(np.vstack(umat), (2,1))
pcaMat = np.tile(np.vstack(pmat),(2,1))

np.savetxt(os.path.join(paths[3], 'verbTopicsManual.txt'), uniqueMat)
np.savetxt(os.path.join(paths[3], 'verbTopicsPCA.txt'), pcaMat)