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


class SKLRegressionMapper(Mapper):
    """NiPy-based GLMMapper implementation

    This is basically a front-end for
    :class:`~ nipy.modalities.fmri.glm.GeneralLinearModel`.
    In particular, it supports all keyword arguments of its
    ``fit()`` method.
    """

    def __init__(self, regs, clf=None, **kwargs):
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
        self.clf = clf


    def _fit_model(self, ds, X, reg_names):
        if self.__clf is None:
            glm = GeneralLinearModel(X)
            glm.fit(ds.samples, **self.glmfit_kwargs)
            out = Dataset(glm.get_beta(),
                          sa={self.get_space(): reg_names})
        else:
            # a model of sklearn linear model
            glm = self.__clf
            glm.fit(X, ds.samples)
            out = Dataset(glm.coef_, sa={self.get_space(): reg_names})
        return glm, out





    def _untrain(self):
        self._transformer = None


    def _get_y(self, ds):
        space = self.get_space()
        if space:
            y = ds.sa[space].value
        else:
            y = None
        return y


    def _get_transformer(self):
        if self._transformer is None:
            self._transformer = deepcopy(self._pristine_transformer)
        return self._transformer


    def _train(self, ds):
        tf = self._get_transformer()
        return tf.fit(ds.samples, self._get_y(ds))


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


    def _forward_data(self, data):
        """Forward-map some data.

        This is a private method that has to be implemented in derived
        classes.

        Parameters
        ----------
        data : anything (supported the derived class)
        """
        raise NotImplementedError


    def _reverse_data(self, data):
        """Reverse-map some data.

        This is a private method that has to be implemented in derived
        classes.

        Parameters
        ----------
        data : anything (supported the derived class)
        """
        raise NotImplementedError