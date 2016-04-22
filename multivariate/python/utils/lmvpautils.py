# utility functions used throughout lmvpa project


def initpaths(platform):
    print "Initializing..."
    import os
    import numpy as np
    p = []
    if 'win' in platform:
        p.append("D:\\fmri\\LanguageMVPA")
        p.append("D:\GitHub\LanguageMVPA\multivariate\python")
    elif 'mac' in platform:
        p.append("/Volumes/fmri/LanguageMVPA")
        p.append('/Users/njchiang/GitHub/LanguageMVPA/multivariate/python')
    elif 'usb' in platform:
        p.append("/Volumes/JEFF/UCLA/LMVPA/")
        p.append('/Users/njchiang/GitHub/LanguageMVPA/multivariate/python')

    p.append(os.path.join(p[1], "labels"))
    p.append(os.path.join(p[0], "Maps", "PyMVPA"))
    # c = ["stim", "verb", "anim", "AP", "CR", "syntax"]

    tmpc = np.loadtxt(os.path.join(p[1], 'labels', 'conds_key.txt'), dtype=str)
    c = {}
    for i, name in enumerate(tmpc[0, :]):
        c[name] = tmpc[1:, i]

    # s = ["LMVPA001", "LMVPA002", "LMVPA003", "LMVPA005", "LMVPA006", "LMVPA007", "LMVPA008", "LMVPA009", "LMVPA010",
    #        "LMVPA011", "LMVPA013", "LMVPA014", "LMVPA015", "LMVPA016", "LMVPA017", "LMVPA018", "LMVPA019"]
    s = {"LMVPA001": ["Run1", "Run2", "Run3", "Run4"],
         "LMVPA002": ["Run1", "Run2", "Run3", "Run4"],
         "LMVPA003": ["Run1", "Run2", "Run3", "Run4", "Run5", "Run6", "Run7", "Run8"],
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

######################################
# load data
def loadsubbetas(p, s, m="grayMatter", a=0):
    # load subject data with paths list, s: sub, c: contrast, m: mask
    print s
    import os
    bsn = str(s + "_LSA_Series.nii.gz")
    bsp = os.path.join(p[0], "betas", "tstat")
    bs = os.path.join(bsp, bsn)
    mnp = os.path.join(p[0], "data", s, "masks")
    mn = str(s+"_"+m+".nii.gz")
    mf = os.path.join(mnp, mn)
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
    if isinstance(c, basestring):
        # must be a list/tuple/array for the logic below
        c = [c]

    # bfn = pjoin(p[0], 'data', s, 'func', 'extra', s+'_'+r+'_mc.nii.gz')
    # motion corrected and coregistered
    bfn = pjoin(p[0], 'data', s, 'func', s + '_' + r + '.nii.gz')
    if m is not None:
        m = pjoin(p[0], 'data', s, 'masks', s+'_'+m+'.nii.gz')
        d = fmri_dataset(bfn, chunks=int(r.split('n')[1]), mask=m)
    else:
        d = fmri_dataset(bfn, chunks=int(r.split('n')[1]), mask=m)
    # This line-- should be different if we're doing GLM, etc.
    efn = pjoin(p[0], 'data', s, 'func', s + '_' + r + '.tsv')
    fe = bids.load_events(efn)
    for ci in c:
        e = adjustevents(fe, ci)
        d = er.assign_conditionlabels(d, e, noinfolabel='rest', label_attr=ci)
    return d


def adjustevents(e, c='trial_type'):
        import numpy as np
        # rounding for now, ONLY because that works for this dataset.
        ee = []
        for i, d in enumerate(e):
            ee.append({'onset': np.round(d['onset']),
                       'duration': np.round(d['duration']),
                       'condition': d[c],
                       'intensity': 1})
        return ee


def replacetargets(d, ckey, c='trial_type'):
    import numpy as np
    if c in ckey:
        d.targets = np.array([ckey[c][np.where(st == ckey['trial_type'])[0][0]] for st in d.sa['trial_type']])
    else:
        print "not a valid contrasts, did not do anything."
    return d


def loadevents(p, s):
    # if isinstance(c, basestring):
    #     # must be a list/tuple/array for the logic below
    #     c = [c]
    fds = {}
    from mvpa2.datasets.sources import bids
    from os.path import join as pjoin
    for sub in s.keys():
        fds[sub] = [bids.load_events(pjoin(p[0], 'data', sub, 'func', sub + '_' + r + '.tsv')) for r in s[sub]]
    return fds


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


def loadmotionparams(p, s):
    import numpy as np
    import os
    res = {}
    for sub in s.keys():
        mcs = [np.loadtxt(os.path.join(p[0], 'data', sub, 'func', 'extra',
                                       sub + '_' + r + '_mc', sub + '_' + r + '_mc.par'))
               for r in s[sub]]
        res[sub] = np.vstack(mcs)
    return res


def amendtimings(ds, b):
    from mvpa2.datasets import eventrelated as er
    import numpy as np
    TR = np.median(np.diff(ds.sa.time_coords))
    idx = 0
    # events are loading wrong...
    theseEvents = b
    events = []
    for i in np.arange(len(b)):
        for ev in theseEvents[i]:
            ev['chunks'] = ds.sa['chunks'].unique[i]
            ev['onset'] += TR*idx
            # ev['targets'] = ev['condition']
            if 'intensity' in ev:
                ev['amplitude'] = ev['intensity']
            else:
                ev['amplitude'] = 1
            if ev['duration'] is not '0':
                events.append(ev)
        idx += np.sum(ds.sa['chunks'].value == ds.sa['chunks'].unique[i])
        # not sure i like the +2, but i think it's right.
        ds.sa.time_coords[ds.sa['chunks'].value == i+2] += ds.sa.time_coords[idx-1] + TR
    return ds, events
#####################################
# regression stuff


def events2dict(events):
    evvars = {}
    for k in events[0]:
        try:
            evvars[k] = [e[k] for e in events]
        except KeyError:
            raise ValueError("Each event property must be present for all "
                             "events (could not find '%s')" % k)
    return evvars

# new make_designmat from scratch
def make_designmat(ds, e, time_attr, condition_attr='targets', design_kwargs=None, regr_attrs=None):
    # make glm regressors for all attributes. so loop through condition_attr and add them all...
    import copy
    from nipy.modalities.fmri.design_matrix import make_dmtx
    import numpy as np
    # Decide/device condition attribute on which GLM will actually be done

    if isinstance(condition_attr, basestring):
        # must be a list/tuple/array for the logic below
        condition_attr = [condition_attr]

    e = copy.deepcopy(e)  # since we are modifying in place
    glm_condition_attrs=[]
    for i, con in enumerate(condition_attr):
        glm_condition_attr = 'regressors_' + str(con)
        glm_condition_attrs.append(glm_condition_attr)
        for ei in e:
            if glm_condition_attr in ei:
                raise ValueError("Event %s already has %s defined.  Should not "
                                 "happen.  Choose another name if defined it"
                                 % (ei, glm_condition_attr))
            ei[glm_condition_attr] = \
                'glm_label_' + str(con) + '+'.join(str(ei[c]) for c in [con])
                # 'glm_label_' + str(con) + '+'.join(str(ei[c]) for c in [con, 'chunks'])

    evvars = events2dict(e)
    add_paradigm_kwargs = {}
    if 'amplitude' in evvars:
        add_paradigm_kwargs['amplitude'] = evvars['amplitude']
    if design_kwargs is None:
        design_kwargs = {}
    if not regr_attrs is None:
        names = []
        regrs = []
        for attr in regr_attrs:
            regr = ds.sa[attr].value
            # add rudimentary dimension for easy hstacking later on
            if regr.ndim < 2:
                regr = regr[:, np.newaxis]
            if regr.shape[1] == 1:
                names.append(attr)
            else:
                #  add one per each column of the regressor
                for i in xrange(regr.shape[1]):
                    names.append("%s.%d" % (attr, i))
            regrs.append(regr)
        regrs = np.hstack(regrs)

        if 'add_regs' in design_kwargs:
            design_kwargs['add_regs'] = np.hstack((design_kwargs['add_regs'],
                                                   regrs))
        else:
            design_kwargs['add_regs'] = regrs
        if 'add_reg_names' in design_kwargs:
            design_kwargs['add_reg_names'].extend(names)
        else:
            design_kwargs['add_reg_names'] = names

    X = {}
    for ci, con in enumerate(condition_attr):
        # create paradigm
        if 'duration' in evvars:
            from nipy.modalities.fmri.experimental_paradigm import BlockParadigm
            # NiPy considers everything with a duration as a block paradigm
            paradigm = BlockParadigm(
                con_id=evvars[glm_condition_attrs[ci]],
                onset=evvars['onset'],
                duration=evvars['duration'],
                **add_paradigm_kwargs)
        else:
            from nipy.modalities.fmri.experimental_paradigm \
                import EventRelatedParadigm
            paradigm = EventRelatedParadigm(
                con_id=evvars[glm_condition_attrs[ci]],
                onset=evvars['onset'],
                **add_paradigm_kwargs)
        X[con] = make_dmtx(ds.sa[time_attr].value, paradigm=paradigm, **design_kwargs)
        for i, reg in enumerate(X[con].names):
            ds.sa[reg] = X[con].matrix[:, i]
        if con in ds.sa.keys():
            ds.sa.pop(con)

        for reg in ds.sa.keys():
            if str(con)+'0' in reg:
                ds.sa['glm_label_probe'] = ds.sa.pop(reg)

    # concatenate X... add chunk regressors...
    # if 'chunks' in ds.sa.keys():
    #     for i in ds.sa['chunks'].unique:
    #         ds.sa['glm_label_chunks' + str(i)] = np.array(ds.sa['chunks'].value == i, dtype=np.int)
    return X, ds

def make_fulldesignmat(ds, e, time_attr, condition_attr='targets', design_kwargs=None, glmfit_kwargs=None, regr_attrs=None):
    # for now... loop through this and feed in one event at a time.
    import copy
    from nipy.modalities.fmri.design_matrix import make_dmtx
    import numpy as np
    # Decide/device condition attribute on which GLM will actually be done
    glm_condition_attr = 'regressor_names'
    if isinstance(condition_attr, basestring):
        # must be a list/tuple/array for the logic below
        condition_attr = [condition_attr]


    glm_condition_attr_map = dict([(con, dict()) for con in condition_attr])  #
    # to map back to original conditions
    e = copy.deepcopy(e)  # since we are modifying in place
    for ei in e:
        if glm_condition_attr in ei:
            raise ValueError("Event %s already has %s defined.  Should not "
                             "happen.  Choose another name if defined it"
                             % (ei, glm_condition_attr))
        compound_label = ei[glm_condition_attr] = \
            'glm_label_' + str(con) + '+'.join(
                str(ei[con]) for con in condition_attr)
        # and mapping back to original values, without str()
        # for each condition:
        for con in condition_attr:
            glm_condition_attr_map[con][compound_label] = ei[con]

    evvars = events2dict(e)
    add_paradigm_kwargs = {}
    if 'amplitude' in evvars:
        add_paradigm_kwargs['amplitude'] = evvars['amplitude']
    # create paradigm
    if 'duration' in evvars:
        from nipy.modalities.fmri.experimental_paradigm import BlockParadigm
        # NiPy considers everything with a duration as a block paradigm
        paradigm = BlockParadigm(
            con_id=evvars[glm_condition_attr],
            onset=evvars['onset'],
            duration=evvars['duration'],
            **add_paradigm_kwargs)
    else:
        from nipy.modalities.fmri.experimental_paradigm \
            import EventRelatedParadigm
        paradigm = EventRelatedParadigm(
            con_id=evvars[glm_condition_attr],
            onset=evvars['onset'],
            **add_paradigm_kwargs)
    # create design matrix -- all kinds of fancy additional regr can be
    # auto-generated
    if design_kwargs is None:
        design_kwargs = {}

    if not regr_attrs is None:
        names = []
        regrs = []
        for attr in regr_attrs:
            regr = ds.sa[attr].value
            # add rudimentary dimension for easy hstacking later on
            if regr.ndim < 2:
                regr = regr[:, np.newaxis]
            if regr.shape[1] == 1:
                names.append(attr)
            else:
                #  add one per each column of the regressor
                for i in xrange(regr.shape[1]):
                    names.append("%s.%d" % (attr, i))
            regrs.append(regr)
        regrs = np.hstack(regrs)

        if 'add_regs' in design_kwargs:
            design_kwargs['add_regs'] = np.hstack((design_kwargs['add_regs'],
                                                   regrs))
        else:
            design_kwargs['add_regs'] = regrs
        if 'add_reg_names' in design_kwargs:
            design_kwargs['add_reg_names'].extend(names)
        else:
            design_kwargs['add_reg_names'] = names

    X = make_dmtx(ds.sa[time_attr].value, paradigm=paradigm, **design_kwargs)
    for i, reg in enumerate(X.names):
        ds.sa[reg] = X.matrix[:, i]

    # for con in condition_attr:
    #     ds.sa.pop(con)
    return X, ds

def corrsig(N, c=None, p=.95):
    # if c exists, this returns the cutoff
    import numpy as np
    from scipy.stats import t
    if not c is None:
        return t.cdf(c/np.sqrt((1-c**2)/(N-2)), N-2)
    else:
        print "functionality not implemented yet, please query a correlation"
        return
        # tsigsq=t.ppf(p, N-2)**2
        # a = tsigsq+N-2
        # return np.sqrt(2*a)/(2*a)

######################################
# classification
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


def testsg(ds, w, p, voxIdx, c='chunks'):
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.signal import savgol_filter
    from mvpa2.mappers.detrend import poly_detrend
    import SavGolFilter as sg
    poly0 = ds.copy(deep=False)
    poly1 = ds.copy(deep=False)
    poly2 = ds.copy(deep=False)
    t = np.arange(ds.shape[0])
    poly_detrend(poly0, polyord=0, chunks_attr=c)
    sgFilt = ds.copy(deep=False)
    manualSet = ds.copy(deep=False)
    manual = getSGDrift(manualSet, w, p)
    # filterMe = filterMat[:,voxIdx]
    # manual = savgol_filter(filterMe, w, p, axis=0)
    # AXIS SHOULD BE 0
    # manualMat = savgol_filter(filterMat, w, p, axis=0)
    poly_detrend(poly1, polyord=1, chunks_attr=c)
    poly_detrend(poly2, polyord=2, chunks_attr=c)
    # sgFilt = savgol_filter(sgFilt, w, p, axis=0)
    sg.sg_filter(sgFilt, window_length=w, polyorder=p, chunks_attr=c, axis=0)
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
            thisreg = sgfilter(d[cinds].copy(deep=False), w, p) - np.mean(d[cinds], axis=0)
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