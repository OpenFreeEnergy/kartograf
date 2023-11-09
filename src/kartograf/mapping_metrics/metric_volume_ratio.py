# This code is part of kartograf and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/kartograf

import logging
import numpy as np
from scipy.spatial import ConvexHull

from gufe.mapping import AtomMapping

from ._abstract_scorer import _AbstractAtomMappingScorer

log = logging.getLogger(__name__)


class MappingVolumeRatioScorer(_AbstractAtomMappingScorer):
    def get_score(self, mapping: AtomMapping) -> float:
        """
            returns a normalized value between 0 and 1, where 0 is the best and 1 ist the worst score.
            The value is rounded to 2 digits.

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
        r = self.get_volume_ratio(mapping)
        return 0.0 if (r < 0) else np.round(r, 2)

    def get_volume_ratio(self, mapping: AtomMapping) -> float:
        """this function calculates the ratio of the volume of the convex hull of the mapped atoms to the volume of the convex hull of the complete molecule

        Parameters
        ----------
        mapping : AtomMapping
            mapped to be scored

        Returns
        -------
        float
            returns the ratio of the volume of the convex hull of the mapped atoms to the volume of the convex hull of the complete molecule

        Raises
        ------
        ValueError
            _description_
        """
        molA = mapping.componentA.to_rdkit()
        molB = mapping.componentB.to_rdkit()
        molA_to_molB = mapping.componentA_to_componentB

        mapping_molA = np.array(list(sorted(molA_to_molB.keys())))
        mapping_molB = np.array(list(sorted(molA_to_molB.values())))

        if len(mapping_molA) < 4 or len(mapping_molB) < 4:
            raise ValueError("Mapping is too small to calculate convex hull")

        complete_molA = ConvexHull(molA.GetConformer().GetPositions()).volume
        map_molA = ConvexHull(
            molA.GetConformer().GetPositions()[mapping_molA]
        ).volume
        complete_molB = ConvexHull(molB.GetConformer().GetPositions()).volume
        map_molB = ConvexHull(
            molB.GetConformer().GetPositions()[mapping_molB]
        ).volume

        ratios = np.array(
            [map_molA / complete_molA, map_molB / complete_molB]
        )
        avg_map_volume_ratio = np.mean(ratios)

        # print("ratios",avg_map_volume_ratio, ratios)
        # print("volumes", map_molA, complete_molA, map_molB, complete_molB)
        # print('ind', mapping_molA, mapping_molB)
        return avg_map_volume_ratio


class MappingRatioMappedAtomsScorer(_AbstractAtomMappingScorer):
    def get_score(self, mapping: AtomMapping) -> float:
        """calculate the number of mapped atoms/number of atoms in the larger molecule

        Parameters
        ----------
        mapping : AtomMapping
            AtomMapping to be scored

        Returns
        -------
        float
            normalized ratio of mapped atoms.
        """
        molA = mapping.componentA.to_rdkit()
        molB = mapping.componentB.to_rdkit()
        molA_to_molB = mapping.componentA_to_componentB

        larger_nAtoms = len(molA.GetAtoms())

        if len(molB.GetAtoms()) > larger_nAtoms:
            larger_nAtoms = len(molB.GetAtoms())

        return np.round((len(molA_to_molB) / larger_nAtoms), 2)
