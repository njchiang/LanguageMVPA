import gensim.models as gm
import os
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

datadir='D:\\fmri\\LanguageMVPA'
labeldir='D:\\GitHub\\LanguageMVPA\\multivariate\\python\\labels'
vecfile='GoogleNews-vectors-negative300.bin.gz'
model = gm.Word2Vec.load_word2vec_format(os.path.join(datadir, vecfile), binary=True)

model.init_sims(replace=True)
words = ['gravy', 'touched', 'broccoli',
         'candle', 'lit', 'match',
         'bus', 'hit', 'truck',
         'tree', 'crushed', 'car',
         'singer', 'kissed', 'guitarist',
         'trainer', 'stretched', 'dancer',
         'goalie', 'kicked', 'referee',
         'woman', 'consoled', 'man']

vecs = []
for i in words:
    vecs.append(model[i])

### get word2vec vectors in shape
a = np.loadtxt(os.path.join(labeldir, 'word2veclabels24.txt'))

sidx = np.arange(0, 24, 3)
subvecs = []
for i in sidx:
    subvecs.append(np.tile(a[i], [4,1]))
vidx = np.arange(1, 24, 3)
verbvecs = []
for i in vidx:
    verbvecs.append(np.tile(a[i], [4,1]))
oidx = np.arange(2, 24, 3)
objvecs = []
for i in oidx:
    objvecs.append(np.tile(a[i], [4,1]))


np.savetxt(os.path.join(labeldir, 'word2vecSubs.txt'), np.repeat(np.vstack(subvecs), 2, axis=0))
np.savetxt(os.path.join(labeldir, 'word2vecVerbs.txt'), np.repeat(np.vstack(verbvecs), 2, axis=0))
np.savetxt(os.path.join(labeldir, 'word2vecObjs.txt'), np.repeat(np.vstack(objvecs), 2, axis=0))
pca = PCA()
pca.fit(np.vstack(vecs))
plt.plot(pca.explained_variance_ratio_)

sentences = ['The gravy touched the broccoli',
             'The candle lit the match',
             'The bus hit the truck',
             'The tree crushed the car',
             'The singer kissed the guitarist',
             'The trainer stretched the dancer',
             'The goalie kicked the referee',
             'The woman consoled the man']
