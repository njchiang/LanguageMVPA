# utility functions used throughout lmvpa project


def initpaths():
    print "Initializing..."
    import os
    p = []
    p.append("D:\\fmri\\LanguageMVPA")
    p.append("D:\GitHub\LanguageMVPA\multivariate\python")
    p.append(os.path.join(p[0], "betas", "tstat"))
    p.append(os.path.join(p[0], "masks", "sub"))
    p.append(os.path.join(p[1], "labels"))
    p.append(os.path.join(p[0], "Maps", "PyMVPA"))
    c = ["stim", "verb", "anim", "AP", "CR", "syntax"]
    # s = ["LMVPA001", "LMVPA002", "LMVPA003", "LMVPA005", "LMVPA006", "LMVPA007", "LMVPA008", "LMVPA009", "LMVPA010",
    #        "LMVPA011", "LMVPA013", "LMVPA014", "LMVPA015", "LMVPA016", "LMVPA017", "LMVPA018", "LMVPA019"]
    s = {"LMVPA001": ["Run1", "Run2", "Run3", "Run4"],
         "LMVPA002": ["Run1", "Run2", "Run3", "Run4"],
         "LMVPA003": ["Run1", "Run2", "Run3", "Run4"],
         "LMVPA005": ["Run1", "Run2", "Run3", "Run4", "Run5", "Run6", "Run7", "Run8"],
         "LMVPA006": ["Run1", "Run2", "Run3", "Run4", "Run5", "Run6", "Run7", "Run8"],
         "LMVPA007": ["Run1", "Run2", "Run3", "Run4", "Run5", "Run6", "Run7", "Run8"],
         "LMVPA008": ["Run1", "Run2", "Run3", "Run4", "Run5", "Run6", "Run7", "Run8"],
         "LMVPA009": ["Run1", "Run2", "Run3", "Run4", "Run5", "Run6", "Run7", "Run8"],
         "LMVPA010": ["Run1", "Run2", "Run3", "Run4", "Run5", "Run6", "Run7", "Run8"],
         "LMVPA011": ["Run1", "Run2", "Run3", "Run4", "Run5", "Run6", "Run7", "Run8"],
         "LMVPA013": ["Run1", "Run2", "Run3", "Run4", "Run5", "Run6", "Run7", "Run8"],
         "LMVPA014": ["Run1", "Run2", "Run3", "Run4", "Run5", "Run6", "Run7", "Run8"],
         "LMVPA015": ["Run1", "Run2", "Run3", "Run4", "Run5", "Run6", "Run7", "Run8"],
         "LMVPA016": ["Run1", "Run2", "Run3", "Run4", "Run5", "Run6", "Run7", "Run8"],
         "LMVPA017": ["Run1", "Run2", "Run3", "Run4", "Run5", "Run6", "Run7", "Run8"],
         "LMVPA018": ["Run1", "Run2", "Run3", "Run4", "Run5", "Run6", "Run7", "Run8"],
         "LMVPA019": ["Run1", "Run2", "Run3", "Run4", "Run5", "Run6", "Run7", "Run8"]
        # "LMVPA020": ["Run1", "Run2", "Run3", "Run4", "Run5", "Run6", "Run7", "Run8"]
           }
    m = ["left_IFG_operc", "left_IFG_triang", "left_STG_post", "left_MTG_post", "grayMatter"]
    return p, s, c, m


def loadsubbetas(p, s, m="grayMatter", a=0):
    # load subject data with paths list, s: sub, c: contrast, m: mask
    print s
    import os
    bsn = str(s + "_LSA_Series.nii.gz")
    bs = os.path.join(p[2], bsn)
    mn = str(s+"_"+m+".nii.gz")
    mf = os.path.join(p[3], mn)
    from mvpa2.datasets.mri import fmri_dataset
    if a != 0:
        return fmri_dataset(samples=bs, targets=a.targets, chunks=a.chunks, mask=mf)
    else:
        return fmri_dataset(samples=bs, mask=mf)


def loadrundata(p, s, r, m=None, c='trial_type'):
    # inputs:
    # p: paths list
    # s: string representing subject ('LMVPA001')
    # r: run ID ('Run1')
    from os.path import join as pjoin
    from mvpa2.datasets import eventrelated as er
    from mvpa2.datasets.mri import fmri_dataset
    from mvpa2.datasets.sources import bids as bids
    efn = pjoin(p[0], 'data', s, 'func', s + '_' + r + '.tsv')
    fe = bids.load_events(efn, 0)
    e = adjustevents(fe, c)
    bfn = pjoin(p[0], 'data', s, 'func', 'extra', s+'_'+r+'_mc.nii.gz')
    if m is not None:
        m = pjoin(p[0], 'data', s, 'masks', s+'_'+m+'.nii.gz')
        d = fmri_dataset(bfn, chunks=int(r.split('n')[1]), mask=m)
    else:
        d = fmri_dataset(bfn, chunks=int(r.split('n')[1]), mask=m)
    return er.assign_conditionlabels(d, e, noinfolabel='rest')


def adjustevents(e, c='trial_type'):
        import numpy as np
        # rounding for now, ONLY because that works for this dataset.
        ee = []
        for i, d in enumerate(e):
            ee.append({'onset': np.round(d['onset']), 'duration': np.round(d['duration']), 'condition': d[c], 'intensity': 1})
        return ee


def loadsubdata(p, s, m=None, c='trial_type'):
    from mvpa2.base import dataset
    fds = {}
    for sub in s.keys():
        print 'loading ' + sub
        rds = [loadrundata(p, sub, r, m, c) for r in s[sub]]
        fds[sub] = dataset.vstack(rds, a=0)
    return fds
# this is PyMVPA's vstack() for merging samples from multiple datasets
# a=0 indicates that the dataset attributes of the first run should be used
# for the merged dataset


def makebinarymodel(l):
    import numpy as np
    m = np.zeros([len(l), len(l)])
    for i in np.unique(l):
        m1 = m.copy()
        m2 = m.copy()
        m1[l == i, :] = 1
        m2[:, l == i] = 1
        m = m1*m2
    return 1-m


def initrsa(c_all, p):
    import numpy as np
    from mvpa2.misc.io.base import SampleAttributes
    import os
    a = np.array([SampleAttributes(os.path.join(p[4], (c + "_attribute_labels.txt"))).targets for c in c_all],
                 dtype=np.int64)
    m = np.zeros([a.shape[1], a.shape[1]])
    for i in np.arange(a.shape[0]):
        m += makebinarymodel(a[i, ])
    from mvpa2.base.learner import ChainLearner
    from mvpa2.mappers.shape import TransposeMapper
    from mvpa2.measures import rsa
    tdsm = rsa.PDistTargetSimilarity(rsa.squareform(m))
    return ChainLearner([tdsm, TransposeMapper()])


######################################
# quality control:

def testsg(ds, w, p, voxIdx):
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.signal import savgol_filter
    from mvpa2.mappers.detrend import poly_detrend
    import SavGolFilter as sg
    poly0 = ds.copy(deep=False)
    poly1 = ds.copy(deep=False)
    poly2 = ds.copy(deep=False)
    t = np.arange(ds.shape[0])
    poly_detrend(poly0, polyord=0, chunks_attr='chunks')
    sgFilt = ds.copy(deep=False)
    manualSet = ds.copy(deep=False)
    manual =  getSGDrift(manualSet, w, p)
    # filterMe = filterMat[:,voxIdx]
    # manual = savgol_filter(filterMe, w, p, axis=0)
    # AXIS SHOULD BE 0
    # manualMat = savgol_filter(filterMat, w, p, axis=0)
    poly_detrend(poly1, polyord=1, chunks_attr='chunks')
    poly_detrend(poly2, polyord=2, chunks_attr='chunks')
    # sgFilt = savgol_filter(sgFilt, w, p, axis=0)
    sg.sg_filter(sgFilt, window_length=w, polyorder=p, chunks_attr='chunks', axis=0)
    plt.plot(t, poly0.samples[:, voxIdx], 'k', label='demean')
    plt.plot(t, poly1.samples[:, voxIdx], 'b+', label='linear')
    plt.plot(t, poly2.samples[:, voxIdx], 'b*', label='quadratic')
    plt.plot(t, sgFilt.samples[:, voxIdx], 'go', label='sg')
    plt.plot(t, manual[:, voxIdx], 'r+', label='sg estimate')
    plt.legend(loc='upper left')
    plt.show()


def getSGDrift(d, w, p, chunks_attr='chunks'):
    import numpy as np
    if chunks_attr is None:
        reg = sgfilter(d.copy(deep=False), w, p) - np.mean(ds, axis=0)
    else:
        uchunks = d.sa[chunks_attr].unique
        reg = []
        for n, chunk in enumerate(uchunks):
            cinds = d.sa[chunks_attr].value == chunk
            thisreg = sgfilter(ds[cinds].copy(deep=False), w, p) - np.mean(ds[cinds], axis=0)
            reg.append(thisreg)
        reg = np.vstack(reg)

    # combine the regs (time x reg)
    # self._regs = np.hstack(reg)
    return np.array(reg)

def sgfilter(data, w, p):
    from scipy.signal import savgol_filter
    mapped = savgol_filter(data.samples, w, p, axis=0)
    return mapped
######################################