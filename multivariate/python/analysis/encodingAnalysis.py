from os.path import join as pjoin
import sys
sys.path.append('D:\\GitHub\\LanguageMVPA\\multivariate\\python\\utils')
sys.path.append('/Users/njchiang/GitHub/LanguageMVPA/multivariate/python/utils')
import lmvpautils as lmvpa
# initialize stuff
paths, subList, contrasts, maskList = lmvpa.initpaths()

# for testing
subjs = {'LMVPA001': subList['LMVPA001']}
ds_all = lmvpa.loadsubdata(paths, subjs, m='left_IFG_operc', c='syntax')
# later this will loop
sub = 'LMVPA001'
ds = ds_all[sub]

from mvpa2.datasets.miscfx import summary
print summary(ds)

import SavGolFilter as sg
b = ds.copy()
sg.sg_filter(b, 79, 3, chunks_attr=None)
b2 = ds.copy()
sg.sg_filter(b2, 79, 3, chunks_attr='chunks')

lmvpa.testsg(ds, 79, 3, 10)

from mvpa2.mappers.detrend import poly_detrend



# savgolfilter=sg.SavGolFilterMapper(window_length=11, polyorder=1, chunks_attr='chunks')


# from mvpa2.mappers.filters import FFTResampleMapper
# resampler=FFTResampleMapper(1, chunks_attr='chunks')
# later will need to discard rest and 0 trials.

# ds.get_mapped(savgolfilter)
# ds.get_mapped(resampler)

# classification analysis below: just trying out.
from mvpa2.mappers.detrend import PolyDetrendMapper
detrender=PolyDetrendMapper(polyord=1, chunks_attr='chunks')
detrended_fds = ds.get_mapped(detrender)

from mvpa2.mappers.zscore import zscore
zscore(detrended_fds, param_est=('targets', ['rest']), chunks_attr='chunks')
fds=detrended_fds

fds = fds[fds.sa.targets != 'rest']
fds = fds[fds.sa.targets != '0']
# from mvpa2.mappers.fx import mean_group_sample
# averager = mean_group_sample(['targets', 'runtype'])
# fds = fds.get_mapped(averager)

from mvpa2.clfs import svm
clf = svm.LinearCSVMC()
from mvpa2.measures.base import CrossValidation
from mvpa2.generators.partition import NFoldPartitioner
from mvpa2.misc import errorfx
cv = CrossValidation(clf,
                     NFoldPartitioner(attr='chunks'),
                     errorfx=errorfx.mean_match_accuracy)

res = cv(fds)