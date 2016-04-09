# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the PyMVPA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
"""Data normalization by Z-Scoring."""

__docformat__ = 'restructuredtext'

import numpy as np

from mvpa2.base import warning
from mvpa2.base.dochelpers import _str, borrowkwargs, _repr_attrs
from mvpa2.mappers.base import accepts_dataset_as_samples, Mapper
from mvpa2.datasets.base import Dataset
from mvpa2.datasets.miscfx import get_nsamples_per_attr, get_samples_by_attr
from mvpa2.support import copy


class EncodingMapper(Mapper):
    """Mapper for running voxel-wise regression with a model"""

    def __init__(self, params=None, param_est=None, chunks_attr='chunks',
                 dtype='float64', **kwargs):
        """
        Parameters
        ----------
        params : None or tuple(mean, std) or dict
          Fixed Z-Scoring parameters (mean, standard deviation). If provided,
          no parameters are estimated from the data. It is possible to specify
          individual parameters for each chunk by passing a dictionary with the
          chunk ids as keys and the parameter tuples as values. If None,
          parameters will be estimated from the training data.
        param_est : None or tuple(attrname, attrvalues)
          Limits the choice of samples used for automatic parameter estimation
          to a specific subset identified by a set of a given sample attribute
          values.  The tuple should have the name of that sample
          attribute as the first element, and a sequence of attribute values
          as the second element. If None, all samples will be used for parameter
          estimation.
        chunks_attr : str or None
          If provided, it specifies the name of a samples attribute in the
          training data, unique values of which will be used to identify chunks of
          samples, and to perform individual Z-scoring within them.
        dtype : Numpy dtype, optional
          Target dtype that is used for upcasting, in case integer data is to be
          Z-scored.
        """

        Mapper.__init__(self, **kwargs)

        self.__chunks_attr = chunks_attr
        self.__params = params
        self.__param_est = param_est
        self.__params_dict = None
        self.__dtype = dtype

        # secret switch to perform in-place z-scoring
        self._secret_inplace = False

    def __repr__(self, prefixes=[]):
        return super(EncodingMapper, self).__repr__(
            prefixes=prefixes
            + _repr_attrs(self, ['params', 'param_est', 'chunks_attr'])
            + _repr_attrs(self, ['dtype'], default='float64'))

