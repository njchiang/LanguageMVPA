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
    c = ["anim", "verb", "syntax", "ActPass", "RelCan", "stimtype", "cross_anim", "cross_verb"]
    s = ["LMVPA001", "LMVPA002", "LMVPA003", "LMVPA005", "LMVPA006", "LMVPA007", "LMVPA008", "LMVPA009", "LMVPA010",
           "LMVPA011", "LMVPA013", "LMVPA014", "LMVPA015", "LMVPA016", "LMVPA017", "LMVPA018", "LMVPA019"]
    m = ["left_IFG_operc", "left_IFG_triang", "left_STG_post", "left_MTG_post", "grayMatter"]
    return p, s, c, m


def loadsub(p, s, m="grayMatter", a=0):
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
