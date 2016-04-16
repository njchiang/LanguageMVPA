# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the PyMVPA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
"""GLMMapper implementation based on the NiPy package."""

__docformat__ = 'restructuredtext'

from mvpa2.base import externals

if externals.exists('nipy', raise_=True):
    from nipy.modalities.fmri.glm import GeneralLinearModel

import numpy as np

from mvpa2.datasets import Dataset
from mvpa2.mappers.glm import GLMMapper


class BootstrapRidgeMapper(GLMMapper):
    """NiPy-based GLMMapper implementation

    This is basically a front-end for
    :class:`~ nipy.modalities.fmri.glm.GeneralLinearModel`.
    In particular, it supports all keyword arguments of its
    ``fit()`` method.
    """

    def __init__(self, regs=[], des=None, glmfit_kwargs=None, **kwargs):
        """
        Parameters
        ----------
        regs : list
          Names of sample attributes to be extracted from an input dataset and
          used as design matrix columns.
        glmfit_kwargs : dict, optional
          Keyword arguments to be passed to GeneralLinearModel.fit().
          By default an AR1 model is used.
        """
        GLMMapper.__init__(self, regs, **kwargs)
        if glmfit_kwargs is None:
            glmfit_kwargs = {}

        self.glmfit_kwargs = glmfit_kwargs

    def _fit_model(self, ds, X, reg_names):
        # X: design matrix. of shape nTimepoints x nRegressors. easy!
        # need to fit per chunk
        glm = BootstrapRidge(X)
        glm.fit(ds.samples, **self.glmfit_kwargs)
        out = Dataset(glm.get_beta(),
                      sa={self.get_space(): reg_names})
        return glm, out


class BootstrapRidge(object):
    """ This class handles the so-called on General Linear Model

    Most of what it does in the fit() and contrast() methods
    fit() performs the standard two-step ('ols' then 'ar1') GLM fitting
    contrast() returns a contrast instance, yileding statistics and p-values.
    The link between fit() and constrast is done vis the two class members:

    glm_results : dictionary of nipy.algorithms.statistics.models.
                 regression.RegressionResults instances,
                 describing results of a GLM fit

    labels : array of shape(n_voxels),
            labels that associate each voxel with a results key
    """

    def __init__(self, X):
        """
        Parameters
        ----------
        X : array of shape (n_time_points, n_regressors)
           the design matrix
        """
        self.X = X
        self.wt_ = None
        self.alpha_ = None
        self.allRcorrs_ = None
        self.valinds_=None

    def fit(self, Y, chunklen=None, nchunks=None, alphas=None, single_alpha=True,
                 normalpha=False, nboots=100, corrmin=.2, singcutoff=1e-10, joined=None,
                 use_corr=True):
        """GLM fitting of a dataset using 'ols' regression or the two-pass

        Parameters
        ----------
        Y : array of shape(n_time_points, n_samples)
            the fMRI data
        model : {'ar1', 'ols'}, optional
            the temporal variance model. Defaults to 'ols'
        steps : int, optional
            Maximum number of discrete steps for the AR(1) coef histogram
        """
        import numpy as np
        if chunklen == None or nchunks == None:
            print 'no chunk length or number of chunks specified... exiting'
            return

        if alphas is None:
            alphas=np.logspace(0,3,20)
        if Y.ndim == 1:
            Y = Y[:, np.newaxis]

        if Y.shape[0] != self.X.shape[0]:
            raise ValueError('Response and predictors are inconsistent')

        # bootstrap ridge regression
        import ridge as ridge
        wt, _, valphas, allRcorrs, valinds  = ridge.bootstrap_ridge(Rstim=self.X,
                                               Rresp=Y,
                                               Pstim=self.X,
                                               Presp=Y,
                                               alphas=alphas,
                                               nboots=nboots,
                                               chunklen=chunklen,
                                               nchunks=nchunks,
                                               corrmin=corrmin,
                                               joined=joined,
                                               singcutoff=singcutoff,
                                               normalpha=normalpha,
                                               single_alpha=single_alpha,
                                               use_corr=use_corr)


        self.wt_ = wt
        self.alpha_ = valphas
        self.allRcorrs_ = allRcorrs
        self.valinds_ = valinds


    def get_beta(self):
        """Accessor for the best linear unbiased estimated of model parameters

        Parameters
        ----------

        Returns
        -------
        beta: array of shape (n_voxels, n_columns)
            the beta
        """
        beta = self.wt_
        return beta

    def get_alpha(self):
        """Accessor for the regularization parameter of the model

        Returns
        -------
        mse: array of shape (n_voxels)
            the sum of square error per voxel
        """
        # build the beta array
        return self.alpha_