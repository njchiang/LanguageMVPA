
def error2acc(d):
    d.samples *= -1
    d.samples += 1
    return d


def run_cv_sl(sl, ds):
    fds = ds.copy(deep=False, sa=['targets', 'chunks'], fa=['voxel_indices'], a=['mapper'])
    # is this necessary if i'm doing tstats? in toolbox it says zscore w.r.t. rest condition,
    # but technically a tstat is already normalized?
    # zscore(fds)
    import time
    wsc_start_time = time.time()
    print "running SL at " + time.strftime("%H:%M:%S")
    thisSL = sl
    res = thisSL(fds)
    print "done in " + str((time.time() - wsc_start_time,)) + " seconds"
    res = error2acc(res)
    return res

# DEPRECATED-- ONLY FOR "SEARCHLIGHTMULTITHREADED"
def run_cc_sl(sl, ds):
    fds = ds.copy(deep=False, sa=['targets', 'chunks'], fa=['voxel_indices'], a=['mapper'])
    # zscore(fds)
    import time
    thisSL = sl
    wsc_start_time = time.time()
    print "running " + str(t) + " at " + time.strftime("%H:%M:%S")
    res = thisSL(fds)
    print "done in " + str((time.time() - wsc_start_time,)) + " seconds"

    res = error2acc(res)
    l2p = res.samples[0]
    p2l = res.samples[1]
    return l2p, p2l
