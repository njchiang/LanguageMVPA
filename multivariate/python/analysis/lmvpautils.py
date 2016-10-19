# utility functions used throughout lmvpa project
zs = lambda v: (v - v.mean(0)) / v.std(0)  ## z-score function


def initpaths(platform):
    print "Initializing..."
    import os
    import numpy as np
    p = []
    if 'win' in platform:
        p.append("D:\\fmri\\LanguageMVPA")
        p.append("D:\\GitHub\\LanguageMVPA\\multivariate\\python")
        p.append("D:\\CloudStation\\Grad\\Research\\LanguageMVPA")
    elif 'mac' in platform:
        p.append("/Volumes/fmri/LanguageMVPA")
        p.append('/Users/njchiang/GitHub/LanguageMVPA/multivariate/python')
        p.append('Users/njchiang/CloudStation/Grad/Research/LanguageMVPA')
    elif 'usb' in platform:
        p.append("/Volumes/JEFF/UCLA/LMVPA/")
        p.append('/Users/njchiang/GitHub/LanguageMVPA/multivariate/python')
        p.append('Users/njchiang/CloudStation/Grad/Research/LanguageMVPA')
    p.append(os.path.join(p[1], "labels"))
    p.append(os.path.join(p[0], "Maps", "PyMVPA"))
    # c = ["stim", "verb", "anim", "AP", "CR", "syntax"]

    tmpc = np.loadtxt(os.path.join(p[1], 'labels', 'conds_key.txt'), dtype=str)
    c = {}
    for i, name in enumerate(tmpc[0, :]):
        c[name] = tmpc[1:, i]

    # manualTopics = np.vstack((np.loadtxt(os.path.join(p[3], 'verbTopicsManual.txt'), dtype=str),
    #                           np.repeat('0', 44), np.repeat('rest', 44)))
    # for i in np.arange(manualTopics.shape[1]):
    #     c['manualTopic' + str(i)] = manualTopics[:,i]
    #
    # pcaTopics = np.vstack((np.loadtxt(os.path.join(p[3], 'verbTopicsPCA.txt'), dtype=str),
    #                        np.repeat('0', 6), np.repeat('rest', 6)))
    # for i in np.arange(pcaTopics.shape[1]):
    #     c['pcaTopic' + str(i)] = pcaTopics[:, i]

    w2vSub = np.vstack((np.loadtxt(os.path.join(p[3], 'word2vecSubs.txt'), dtype=str),
                             np.repeat('0', 300), np.repeat('rest', 300)))
    w2vVerb = np.vstack((np.loadtxt(os.path.join(p[3], 'word2vecVerbs.txt'), dtype=str),
                             np.repeat('0', 300), np.repeat('rest', 300)))
    w2vObj = np.vstack((np.loadtxt(os.path.join(p[3], 'word2vecObjs.txt'), dtype=str),
                             np.repeat('0', 300), np.repeat('rest', 300)))
    np.random.seed(100)
    randfeat = np.vstack((np.random.rand(64, 302).astype(str),
                         np.repeat('0', 302), np.repeat('rest', 302)))

    w2vFeatures = w2vVerb
    for i in np.arange(w2vFeatures.shape[1]):
        c['word2vec' + str(i)] = w2vFeatures[:, i]

    for i in np.arange(randfeat.shape[1]):
        c['random' + str(i)] = randfeat[:, i]
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
def loadsubbetas(p, s, btype='tstat', m="grayMatter"):
    # load subject data with paths list, s: sub, c: contrast, m: mask
    print s
    import os
    bsp = os.path.join(p[0], "data", s, "func")
    bsn = str(s + "_LSS_" + btype + ".nii.gz")
    bs = os.path.join(bsp, bsn)
    mnp = os.path.join(p[0], "data", s, "masks")
    mn = str(s + "_" + m + ".nii.gz")
    mf = os.path.join(mnp, mn)
    from mvpa2.datasets.mri import fmri_dataset
    fds = fmri_dataset(samples=bs, mask=mf)
    import pandas as pd
    attrs = pd.read_csv(os.path.join(bsp, str(s+"_LSS_betas.tsv")), sep='\t')
    fds.sa['chunks'] = attrs['run'].tolist()
    for c in attrs.keys():
        fds.sa[c] = attrs[c].tolist()
    return fds


def loadrundata(p, s, r, m=None, c=None):
    # inputs:
    # p: paths list
    # s: string representing subject ('LMVPA001')
    # r: run ID ('Run1')
    from os.path import join as pjoin
    from mvpa2.datasets import eventrelated as er
    from mvpa2.datasets.mri import fmri_dataset
    from mvpa2.datasets.sources import bids as bids


    # bfn = pjoin(p[0], 'data', s, 'func', 'extra', s+'_'+r+'_mc.nii.gz')
    # motion corrected and coregistered
    bfn = pjoin(p[0], 'data', s, 'func', s + '_' + r + '.nii.gz')
    if m is not None:
        m = pjoin(p[0], 'data', s, 'masks', s+'_'+m+'.nii.gz')
        d = fmri_dataset(bfn, chunks=int(r.split('n')[1]), mask=m)
    else:
        d = fmri_dataset(bfn, chunks=int(r.split('n')[1]))
    # This line-- should be different if we're doing GLM, etc.
    efn = pjoin(p[0], 'data', s, 'func', s + '_' + r + '.tsv')
    fe = bids.load_events(efn)
    if c is None:
        tmpe = events2dict(fe)
        c = tmpe.keys()
    if isinstance(c, basestring):
        # must be a list/tuple/array for the logic below
        c = [c]
    for ci in c:
        e = adjustevents(fe, ci)
        d = er.assign_conditionlabels(d, e, noinfolabel='rest', label_attr=ci)
    return d


def adjustevents(e, c='trial_type'):
        import numpy as np
        # rounding for now, ONLY because that works for this dataset. But should not round for probe
        ee = []
        for i, d in enumerate(e):
            if int(d['probe']) == 0:
                ee.append({'onset': np.round(d['onset']),
                           'duration': np.round(d['duration']),
                           'condition': d[c],
                           'intensity': int(d['intensity'])})
            else:
                ee.append({'onset': d['onset'],
                           'duration': d['duration'],
                           'condition': d[c],
                           'intensity': int(d['intensity'])})
        return ee


def replacetargets(d, ckey, c='trial_type'):
    import numpy as np
    if c in ckey:
        d.sa[c] = [ckey[c][np.where(st == ckey['trial_type'])[0][0]] for st in d.sa['trial_type']]
        d.sa['targets'] = d.sa[c]
    else:
        print "not a valid contrasts, did not do anything."
    return d


def sortds(d, c='trial_type'):
    ds = d.copy()
    idx = [i[0] for i in sorted(enumerate(ds.sa[c].value), key=lambda x: x[1])]
    ds = d[idx, :]
    return ds


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


def loadsubdata(p, s, m=None, c=None):
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


def amendtimings(ds, b, extras=None):
    # have this add all of the extra stuff
    from mvpa2.datasets import eventrelated as er
    import numpy as np
    TR = np.median(np.diff(ds.sa.time_coords))
    idx = 0
    # events are loading wrong...
    theseEvents = b
    events = []
    for i, te in enumerate(theseEvents):
        for ev in te:
            ev['chunks'] = ds.sa['chunks'].unique[i]
            ev['onset'] += idx
            # ev['targets'] = ev['condition']
            if 'intensity' in ev:
                ev['amplitude'] = ev['intensity']
            else:
                ev['amplitude'] = 1
            # add extra regressors
            if extras is not None:
                for k in extras.keys():
                    if ('manualTopic' in k) or ('pcaTopic' in k) or \
                            ('word2vec' in k) or ('random' in k):
                        ev[k] = extras[k][extras['trial_type'] == ev['trial_type']][0]
            if ev['duration'] is not '0':
                events.append(ev)

        if i < len(b)-1:
            # DO I NEED TO ADD TR TO THE FIRST ONE?!?!
            ds.sa['time_coords'].value[ds.chunks == i+2] += np.max(ds.sa['time_coords'][ds.chunks == i+1]) + TR
            idx = np.min(ds.sa['time_coords'][ds.chunks == i+2])

    # NOT SURE IF I SHOULD DO THIS OR NOT... NO
    # ds.sa['time_coords'].value += TR
    return ds, events


def events2dict(events):
    evvars = {}
    for k in events[0]:
        try:
            evvars[k] = [e[k] for e in events]
        except KeyError:
            raise ValueError("Each event property must be present for all "
                             "events (could not find '%s')" % k)
    return evvars


#####################################
# regression stuff
def condensematrix(dm, pd, names, key, hrf='canonical', op='mult'):
    # returns condition with probe removed
    import copy as cp
    import numpy as np
    delays = None
    if hrf == 'fir':
        delays = []
        for i in dm.names:
            if i == 'constant':
                delays.append('-1')
            else:
                delays.append(i.split('_')[i.split('_').index('delay') + 1])
        delays = np.array(delays, dtype=int)

    if op == 'stack':
        for i in dm.names:
            if (i != 'constant'):
                if (i.split('_')[i.split('_').index(key) + 1] != '0'):
                    pd.append(dm.matrix[:, dm.names.index(i)])
                    names.append(i.replace('glm_label_', ''))
    else:
        idx = []
        for i in dm.names:
            if i == 'constant':
                idx.append('0')
            else:
                idx.append(i.split('_')[i.split('_').index(key)+1])
        idx = np.array(idx, dtype=float)

        if delays is not None:
            for d in np.arange(np.max(delays)+1):
                outkey = key + '_delay_' + str(d)
                outidx = idx[delays == d]
                pd.append(np.dot(dm.matrix[:, delays==d], outidx))
                names.append(outkey)
        else:
            pd.append(np.dot(dm.matrix, idx))
            names.append(key)

# new make_designmat from scratch
def make_designmat(ds, eorig, time_attr='time_coords', condition_attr='targets', design_kwargs=None, regr_attrs=None):
    # make glm regressors for all attributes. so loop through condition_attr and add them all...
    import copy
    from nipy.modalities.fmri.design_matrix import make_dmtx
    import numpy as np
    # Decide/device condition attribute on which GLM will actually be done
    if isinstance(condition_attr, basestring):
        # must be a list/tuple/array for the logic below
        condition_attr = [condition_attr]

    e = copy.deepcopy(eorig)  # since we are modifying in place
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
                'glm_label_' + str(con) + '_' + '+'.join(str(ei[c]) for c in [con])

    evvars = events2dict(e)
    add_paradigm_kwargs = {}
    if 'amplitude' in evvars:
        add_paradigm_kwargs['amplitude'] = evvars['amplitude']
    if design_kwargs is None:
        design_kwargs = {}
    if regr_attrs is not None:
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


def make_parammat(dm, hrf='canonical', zscore=False):
    # remove anything with a 0 and include probe as a feature
    # assuming dm is a dict
    import numpy as np
    out = dm[dm.keys()[0]]
    pd = []
    names = []
    for key in dm.keys():
        if key == 'motion':
            names.append('motion_0')
            pd.append(np.dot(dm[key].matrix, np.array([1, 0, 0, 0, 0, 0, 0])))
            names.append('motion_1')
            pd.append(np.dot(dm[key].matrix, np.array([0, 1, 0, 0, 0, 0, 0])))
            names.append('motion_2')
            pd.append(np.dot(dm[key].matrix, np.array([0, 0, 1, 0, 0, 0, 0])))
            names.append('motion_3')
            pd.append(np.dot(dm[key].matrix, np.array([0, 0, 0, 1, 0, 0, 0])))
            names.append('motion_4')
            pd.append(np.dot(dm[key].matrix, np.array([0, 0, 0, 0, 1, 0, 0])))
            names.append('motion_5')
            pd.append(np.dot(dm[key].matrix, np.array([0, 0, 0, 0, 0, 1, 0])))
        # hardcode stim and verb
        elif key == 'stim' or key == 'verb' or key == 'anim':
            condensematrix(dm[key], pd, names, key, hrf, op='stack')
        else:
            condensematrix(dm[key], pd, names, key, hrf, op='mult')
    # don't need constant because normalized data
    # pd.append(np.ones(np.shape(pd[-1])))
    # names.append('constant')
    if zscore==True:
        out.matrix = zs(np.array(pd).T)
    else:
        out.matrix = (np.array(pd).T)
    out.names = names
    return out


def testmodel(wts, des, ds, tc, use_corr=True):
    import numpy as np
    widx = wts.sa['chunks'].unique
    didx = ds.sa['chunks'].unique
    if len(widx) != len(didx):
        print "unequal number of chunks... exiting"
        return
    if 'word2vec' in tc:
        tc.remove('word2vec')
        for i in np.arange(0, 300):
            tc.append('word2vec' + str(i))

    corrs=[]
    regidx = [des.names.index(i) for i in tc]
    for i in np.arange(len(widx)):
        pred = np.dot(des.matrix[:, regidx],
                      wts[wts.sa['chunks'].value == widx[i]].samples[regidx, :])[ds.sa['chunks'].value == didx[i]]
        Presp = ds[ds.sa['chunks'].value == didx[i]].samples
        # Find prediction correlations
        nnpred = np.nan_to_num(pred)
        if use_corr:
            vcorrs = np.nan_to_num(np.array([np.corrcoef(Presp[:, ii], nnpred[:, ii].ravel())[0, 1]
                                            for ii in range(Presp.shape[1])]))
        else:
            resvar = (Presp - pred).var(0)
            Rsqs = 1 - (resvar / Presp.var(0))
            vcorrs = np.sqrt(np.abs(Rsqs)) * np.sign(Rsqs)
        corrs.append(vcorrs)
    from mvpa2.datasets import Dataset
    return Dataset(np.vstack(corrs), sa={'chunks': ds.sa['chunks'].unique}, fa=ds.fa, a=ds.a)


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