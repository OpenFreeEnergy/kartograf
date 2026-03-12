# This code is part of kartograf and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/kartograf

import abc
import logging

from gufe.mapping import AtomMapping

logger = logging.getLogger(__name__)


class _AbstractAtomMappingScorer(abc.ABC):
    def __init__(self) -> None:
        pass

    def __call__(self, mapping: AtomMapping, *args, **kwargs) -> float:
        """Call the scorer to get the score of an atom mapping.

        Parameters
        ----------
        mapping : AtomMapping
            The atom mapping to be scored.
        *args
            Additional positional arguments.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        float
            A score between 0 and 1, where 0 indicates a very good score and 1 a very bad score.
        """
        return self.get_score(mapping)

    @abc.abstractmethod
    def get_score(self, mapping: AtomMapping, *args, **kwargs) -> float:
        """Calculate the score of an atom mapping.

        The scoring function returns a value between 0 and 1, where a value close to 1.0 indicates a small distance,
        and a score close to zero indicates a large cost/error.

        Parameters
        ----------
        mapping : AtomMapping
            The atom mapping to be scored.
        *args
            Additional positional arguments.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        float
            A score between 0 and 1, where 0 indicates a very good score and 1 a very bad score.
        """
        pass
