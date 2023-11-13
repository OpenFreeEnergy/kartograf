# This code is part of kartograf and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/kartograf

import logging
import numpy as np
from scipy import constants as const

from gufe.mapping import AtomMapping

from ._abstract_scorer import _AbstractAtomMappingScorer

logger = logging.getLogger(__name__)


class MappingRMSDScorer(_AbstractAtomMappingScorer):
    def get_rmsd(self, mapping: AtomMapping) -> float:
        """this function calculates the rmsd between
        the mapped atoms of the two molecules

        Parameters
        ----------
        mapping : AtomMapping
            A mapping to be scored

        Returns
        -------
        float
            returns non-normalized rmsd
        """
        molA = mapping.componentA.to_rdkit()
        molB = mapping.componentB.to_rdkit()
        molA_to_molB = mapping.componentA_to_componentB

        molPosA = molA.GetConformer().GetPositions()
        molPosB = molB.GetConformer().GetPositions()

        diffs = []
        for i, j in molA_to_molB.items():
            diff = molPosB[j] - molPosA[i]
            diffs.append(diff)

        diffs = np.array(diffs)
        rmsd_map_diff = np.sqrt(np.sum(diffs**2))

        return float(rmsd_map_diff)

    def get_rmsd_p(
        self,
        mapping: AtomMapping,
        accepted_distance_rmsd: float = 0.5,
        k_hook: float = 1,
        T: float = 298,
    ) -> float:
        """estimate likelihood of this shift by calculating the probability
        of the rmsd of the mapping with a harmonic oscillator

        Parameters
        ----------
        mapping : AtomMapping
            A mapping to be scored
        k_hook : float, optional
            hook constant of the harmonic oscillator
        T : float, optional
            Temperature
        Returns
        -------
        float
            likelihood of the shift (1-p)
        """
        rmsd = self.get_rmsd(mapping)
        beta = 1000 / (const.k * const.Avogadro * T)
        V = k_hook * (rmsd - accepted_distance_rmsd)
        p = np.exp(-beta * V) if (np.exp(-beta * V) < 1) else 1
        return float(p)

    def get_score(
        self,
        mapping,
        accepted_distance_rmsd: float = 0,
        k_hook: float = 1,
        T: float = 298,
    ) -> float:
        """ Calculate mapping RMSD based score
            returns a normalized value between 0 and 1, where 1.0 is the best
            and 0.0 is the worst score.

        Parameters
        ----------
        mapping : AtomMapping
            A mapping to be scored
        k_hook : float, optional
            hook constant of the harmonic oscillator
        T : float, optional
            Temperature
        Returns
        -------
        float
            normalized score between 0 and 1.
        """
        s = self.get_rmsd_p(mapping,
                            accepted_distance_rmsd=accepted_distance_rmsd,
                            k_hook=k_hook, T=T,)

        return float(s)
