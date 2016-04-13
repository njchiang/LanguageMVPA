# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the PyMVPA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
"""Transform consecutive samples into individual multi-dimensional samples"""

__docformat__ = 'restructuredtext'

import numpy as np

from mvpa2.mappers.base import Mapper
from mvpa2.clfs.base import accepts_dataset_as_samples
from mvpa2.base.dochelpers import _str
import copy
if __debug__:
    from mvpa2.base import debug

def _assure_consistent_a(ds, oshape):
    """If ds.shape differs from oshape, invoke set_length_check
       for the corresponding collection
    """
    shape = ds.shape
    if oshape[0] != shape[0]:
        ds.sa.set_length_check(shape[0])
    if oshape[1] != shape[1]:
        ds.fa.set_length_check(shape[1])


class BootstrapRidge(Mapper):
    """Mapper to combine multiple samples into a single sample.

    Notes
    -----

    This mapper is somewhat unconventional since it doesn't preserve number
    of samples (ie the size of 0-th dimension).
    """
    # TODO: extend with the possibility to provide real onset vectors and a
    #       samples attribute that is used to determine the actual sample that
    #       is matching a particular onset. The difference between target onset
    #       and sample could be stored as an additional sample attribute. Some
    #       utility functionality (outside BoxcarMapper) could be used to merge
    #       arbitrary sample attributes into the samples matrix (with
    #       appropriate mapper adjustment, e.g. CombinedMapper).
    def __init__(self, des=None, alphas=None, chunklen=None, nchunks=None, single_alpha=True,
                 normalpha=False, nboots=100, corrmin=.2, singcutoff=1e-10, joined=None,
                 use_corr=True, **kwargs):
        """
        Parameters
        ----------
        startpoints : sequence
          Index values along the first axis of 'data'.
        boxlength : int
          The number of elements after 'startpoint' along the first axis of
          'data' to be considered for the boxcar.
        offset : int
          The offset between the provided starting point and the actual start
          of the boxcar.
        """
        Mapper.__init__(self, **kwargs)
        self._des = des
        if alphas is None:
            self._alphas = np.logspace(0,3,20)
        else:
            self._alphas = alphas
        self._sa_filter = 'chunks'
        self._fa_filter = None
        self._a_filter = None
        self._nboots=nboots
        self._chunklen=chunklen
        self._nchunks=nchunks
        self._corrmin=corrmin
        self._joined=joined
        self._singcutoff=singcutoff
        self._normalpha = normalpha
        self._single_alpha = single_alpha
        self._use_corr = use_corr


    def _forward_dataset(self, dataset):
        """Forward-map a dataset.

        This is a private method that can be reimplemented in derived
        classes. The default implementation forward-maps the dataset samples
        and returns a new dataset that is a shallow copy of the input with
        the mapped samples.

        Parameters
        ----------
        dataset : Dataset-like
        """
        if __debug__:
            debug('MAP_', "Forward-map %s-shaped samples in dataset with '%s'."
                        % (dataset.samples.shape, self))
        msamples, wts, alphas = self._forward_data(dataset.samples, dataset.sa['chunks'])
        if __debug__:
            debug('MAP_', "Make shallow copy of to-be-forward-mapped dataset "
                    "and assigned forward-mapped samples ({sf}a_filters: "
                    "%s, %s, %s)." % (self._sa_filter, self._fa_filter,
                                      self._a_filter))
        mds = dataset.copy(deep=False,
                           sa=self._sa_filter,
                           fa=self._fa_filter,
                           a=self._a_filter)

        mds.samples = msamples
        mds.sa.wts = wts
        mds.sa.alphas = alphas
        # _assure_consistent_a(mds, dataset.shape)

        if __debug__:
            debug('MAP_', "Return forward-mapped dataset.")
        return mds

    def _forward_data(self, data, idx):
        des = self._des
        samples = []
        # _assure_consistent_a(mds, dataset.shape)
        import ridge as ridge
        wts = []
        alphas = []
        for c in np.arange(len(idx.unique)):
            testidx = idx.value==idx.unique[c]
            trainidx = idx.value!=idx.unique[c]
            wt, corrs, valphas, allRcorrs, valinds = ridge.bootstrap_ridge(Rstim=des.matrix[trainidx],
                                                                   Rresp=data[trainidx],
                                                                   Pstim=des.matrix[testidx],
                                                                   Presp=data[testidx],
                                                                   alphas=self._alphas,
                                                                   nboots=self._nboots,
                                                                   chunklen=self._chunklen,
                                                                   nchunks=self._nchunks,
                                                                   corrmin=self._corrmin,
                                                                   joined=self._joined,
                                                                   singcutoff=self._singcutoff,
                                                                   normalpha=self._normalpha,
                                                                   single_alpha=self._single_alpha,
                                                                   use_corr=self._use_corr)

            samples.append(corrs)
            wts.append(wt)
            alphas.append(valphas)
        return samples, wts, alphas

def bootstrap_ridge(ds, *args, **kwargs):
    """IIR-based frequency filtering.

    Parameters
    ----------
    ds : Dataset
    **kwargs
      For all other arguments, please see the documentation of
      IIRFilterMapper.
    """
    dm = BootstrapRidge(*args, **kwargs)
    # dm._secret_inplace_filter = True
    # map
    return dm.forward(ds)
    # and append the mapper to the dataset
    # mapped._append_mapper(dm)
    #  for now it's not appended... but just keep that in mind.

