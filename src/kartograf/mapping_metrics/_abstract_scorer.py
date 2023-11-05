# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/kartograf

import abc
import logging
import numpy as np

from gufe.mapping import AtomMapping

logger = logging.getLogger(__name__)

eukli = lambda x, y: np.sqrt(np.sum(np.square(y - x)))
rms_func = lambda x: np.sqrt(np.mean(np.square(x)))


class _AbstractAtomMappingScorer(abc.ABC):
    def __init__(self):
        pass

    def __call__(self, mapping: AtomMapping, *args, **kwargs) -> float:
        return self.get_score(mapping)

    @abc.abstractmethod
    def get_score(self, mapping: AtomMapping, *args, **kwargs) -> float:
        """
            the scoring function returns a value between 0 and 1.
            a value close to 1.0 indicates a small distance, a score close to zero indicates a large cost/error.

        Parameters
        ----------
        mapping: AtomMapping
            the mapping to be scored
        args
        kwargs

        Returns
        -------
        float
            a value between [0,1] where one is a very bad score and 0 a very good one.

        """
        pass
