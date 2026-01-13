# This code is part of kartograf and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/kartograf

import abc
import logging

from gufe.mapping import AtomMapping

logger = logging.getLogger(__name__)


class _AbstractAtomMappingScorer(abc.ABC):
    def __init__(self):
        pass

    def __call__(self, mapping: AtomMapping, *args, **kwargs) -> float:
        return self.get_score(mapping)

    @abc.abstractmethod
    def get_score(self, mapping: AtomMapping, *args, **kwargs) -> float:
        """calculate the score
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
