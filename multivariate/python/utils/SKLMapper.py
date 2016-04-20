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
import numpy as np
from mvpa2.support.copy import deepcopy
from mvpa2.datasets import Dataset
from mvpa2.mappers.base import Mapper


class SKLRegressionMapper(Mapper):
    """NiPy-based GLMMapper implementation

    This is basically a front-end for
    :class:`~ nipy.modalities.fmri.glm.GeneralLinearModel`.
    In particular, it supports all keyword arguments of its
    ``fit()`` method.
    """

    def __init__(self, regs=[], add_regs=None, clf=None, **kwargs):
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
        Mapper.__init__(self, auto_train=True, **kwargs)
        # GLMMapper.__init__(self, regs, **kwargs)
        self._clf = None
        self._pristine_clf = clf
        self.regs = list(regs)
        if add_regs is None:
            add_regs = tuple()
        self.add_regs = tuple(add_regs)

    def _build_design(self, ds):
        X = None
        regsfromds = list(self.regs)
        reg_names=None
        if len(regsfromds):
            X = np.vstack([ds.sa[reg].value for reg in regsfromds]).T
            reg_names=regsfromds
        if len(self.add_regs):
            regs=[]
            if reg_names is None:
                reg_names = []
            for reg in self.add_regs:
                regs.append(reg[1])
                reg_names.append(reg[0])
            if X is None:
                X = np.vstack(regs).T
            else:
                X = np.vstack([X.T] + regs).T
        if self.params.add_constant:
            constant = np.ones(len(ds))
            if X is None:
                X = constant[None].T
            else:
                X = np.vstack((X.T, constant)).T
            if reg_names is None:
                reg_names = ['constant']
            else:
                reg_names.append('constant')
        if X is None:
            raise ValueError("no design specified")
        return reg_names, X

    def _fit_model(self, ds, X, reg_names):
        # a model of sklearn linear model
        glm = self.__clf
        glm.fit(X, ds.samples)
        out = Dataset(glm.coef_, sa={self.get_space(): reg_names})
        return glm, out

    def _get_y(self, ds):
        space = self.get_space()
        if space:
            y = ds.sa[space].value
        else:
            y = None
        return y

    def _get_clf(self):
        if self._clf is None:
            self._clf = deepcopy(self._clf)
        return self._clf

    # def _train(self, ds):
    #     tf = self._get_clf()
    #     reg_names, X = self._build_design(self, ds)
    #     return tf.fit(X, ds.samples)


    # def _forward_dataset(self, ds):
    #     tf = self._get_transformer()
    #     if not self.is_trained:
    #         # sklearn support fit and transform at the same time, which might
    #         # be a lot faster, but we only do that, if the mapper is not
    #         # trained already
    #         out = tf.fit_transform(ds.samples, self._get_y(ds))
    #         self._set_trained()
    #     else:
    #         # some SKL classes do not swallow a superfluous `y` argument
    #         # we could be clever and use 'inspect' to query the function
    #         # signature, but we'll use a sledge hammer instead
    #         try:
    #             out = tf.transform(ds.samples, self._get_y(ds))
    #         except TypeError:
    #             out = tf.transform(ds.samples)
    #     return out

    # I don't want multivariate regression... do I? I want iterated regression... okay so it works. output is feature x beta

    def _forward_data(self, data):
        """Forward-map some data instead of implementing forward_dataset (which will call this on a copy).

        This is a private method that has to be implemented in derived
        classes.

        Parameters
        ----------
        data : anything (supported the derived class)
        """
        reg_names, X = self._build_design(self, data)
        model, out = self._fit_model(data, X, reg_names)
        out.fa.update(data.fa)
        out.a.update(data.a) # this last one might be a bit to opportunistic
        # determine the output
        if self.params.return_design:
            if not len(out) == len(X.T):
                raise ValueError("cannot include GLM regressors as sample "
                                 "attributes (dataset probably contains "
                                 "something other than parameter estimates")
            out.sa['regressors'] = X.T
        if self.params.return_model:
            out.a['model'] = model


    def _reverse_data(self, data):
        """Reverse-map some data.

        This is a private method that has to be implemented in derived
        classes.

        Parameters
        ----------
        data : anything (supported the derived class)
        """
        raise NotImplementedError