# prototype...
import itertools
import numpy as np
from scipy import stats
import pylab as pl
from sklearn import svm, linear_model, cross_validation

def ranktransform(x, y):
# apply ranking transform
# form all pairwise combinations
    comb = itertools.combinations(range(x.shape[0]), 2)
    k = 0
    Xp, yp, diff = [], [], []
    for (i, j) in comb:
        if y[i] == y[j]:
            # skip if same target or different group
            continue
        Xp.append(x[i] - x[j])
        diff.append(y[i] - y[j])
        yp.append(np.sign(diff[-1]))
        # output balanced classes
        if yp[-1] != (-1) ** k:
            yp[-1] *= -1
            Xp[-1] *= -1
            diff[-1] *= -1
        k += 1
    Xp, yp, diff = map(np.asanyarray, (Xp, yp, diff))

    return Xp, yp, diff



clf = svm.SVC(kernel='linear', C=.1)
clf.fit(Xp, yp)
coef = clf.coef_.ravel() / linalg.norm(clf.coef_)
