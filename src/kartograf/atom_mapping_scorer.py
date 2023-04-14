# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/kartograf

import math
import numpy as np
from typing import Tuple

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdShapeHelpers

from scipy.spatial import ConvexHull
from scipy import constants as const

from gufe import LigandAtomMapping
from gufe.mapping import AtomMapping

import logging

log = logging.getLogger(__name__)

"""
    Metrics
"""
# simple metrics
eukli = lambda x, y: np.sqrt(np.sum(np.square(y - x)))
rms_func = lambda x: np.sqrt(np.mean(np.square(x)))


import abc
class _abstract_ligand_atom_mapping_scorer(abc.ABC):

    def __call__(self, mapping:LigandAtomMapping, *args, **kwargs)->float:
        return self.get_score(mapping)

    @abc.abstractmethod
    def get_score(self,mapping:LigandAtomMapping, *args, **kwargs)->float:
        """
            the scoring function returns a value between 0 and 1.
            a value close to Zero indicates a small distance, a score close to one indicates a large cost/error.

        Parameters
        ----------
        mapping: LigandAtomMapping
            the mapping to be scored
        args
        kwargs

        Returns
        -------
        float
            a value between [0,1] where one is a very bad score and 0 a very good one.

        """
        pass

class mapping_rmsd(_abstract_ligand_atom_mapping_scorer):
    def get_rmsd(self, mapping: LigandAtomMapping) -> float:
        """this function calculates the rmsd between the mapped atoms of the two molecules

        Parameters
        ----------
        mapping : LigandAtomMapping
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
        rmsd_map_diff = np.round(np.sqrt(np.sum(diffs**2)), 3)

        return np.round(rmsd_map_diff,5)

    def get_rmsd_p(self, mapping, accepted_distance_rmsd:float = 0.5, k_hook=1, T=298) -> float:
        """estimate likelihood of this shift by calculating the probability of the rmsd of the mapping with a harmonic oscillator

        Parameters
        ----------
        mapping : _type_
            _description_
        k_hook : float, optional
            hook constant of the harmonic oscillator
        T : float, optional
            Temperature
        Returns
        -------
        float
            1 - likelihood of the shift (1-p)
        """
        rmsd = self.get_rmsd(mapping)
        beta = 1000 / (const.k * const.Avogadro * T)
        V = k_hook * (rmsd-accepted_distance_rmsd)
        p = np.exp(-beta * V) if(np.exp(-beta * V) < 1) else 1
        return np.round(p,5)

    def get_score(self, mapping, accepted_distance_rmsd:float = 0, k_hook=10**0, T=298)->float:
        return np.round(1-self.get_rmsd_p(mapping, accepted_distance_rmsd = accepted_distance_rmsd, k_hook= k_hook, T= T),2)


class mapping_volume_ratio(_abstract_ligand_atom_mapping_scorer):

    def get_score(self, mapping: LigandAtomMapping) -> float:
        return np.round(1-self.get_volume_ratio(mapping), 2)

    def get_volume_ratio(self, mapping: LigandAtomMapping) -> float:
        """this function calculates the ratio of the volume of the convex hull of the mapped atoms to the volume of the convex hull of the complete molecule

        Parameters
        ----------
        mapping : LigandAtomMapping
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

        if len(mapping_molA) < 4:
            return np.inf
            raise ValueError("Mapping is too small to calculate convex hull")

        complete_molA = ConvexHull(molA.GetConformer().GetPositions()).volume
        map_molA = ConvexHull(molA.GetConformer().GetPositions()[mapping_molA]).volume
        complete_molB = ConvexHull(molB.GetConformer().GetPositions()).volume
        map_molB = ConvexHull(molB.GetConformer().GetPositions()[mapping_molB]).volume

        ratios = np.array([map_molA / complete_molA, map_molB / complete_molB])
        avg_map_volume_ratio = np.mean(ratios)

        print("ratios",avg_map_volume_ratio, ratios)
        print("volumes", map_molA, complete_molA, map_molB, complete_molB)
        #print('ind', mapping_molA, mapping_molB)
        return avg_map_volume_ratio


class mapping_ratio_of_mapped_atoms(_abstract_ligand_atom_mapping_scorer):
    def get_score(self, mapping: LigandAtomMapping) -> float:
        """calculate the number of mapped atoms/number of atoms in the larger molecule

        Parameters
        ----------
        mapping : LigandAtomMapping
            ligandAtomMapping to be scored

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

        return np.round(1-(len(molA_to_molB) / larger_nAtoms), 2)


class mapping_shape_distance(_abstract_ligand_atom_mapping_scorer):
    mapping_mols: Tuple[Chem.Mol, Chem.Mol]

    def get_score(self, mapping:LigandAtomMapping) ->float:
        s = self.get_shape_mapping_distance(mapping)
        return np.round(s,2) if(s<1) else 1

    def get_shape_mapping_distance(self, mapping:LigandAtomMapping, _rd_shape_dist_func=rdShapeHelpers.ShapeTanimotoDist,
                                   _gridSpacing=0.5, _vdwScale=0.8, _ignoreHs=False, _maxLayers=-1,
                                   _stepSize=0.25)->float:
        """
            the function builds two mapped mols, each for one molecule. Basically reduces the atoms to the mapped region.
            next it calculates the shape distance of both full molecules and of both mapped mols.
            The ratio mappedMolShapeDist/MolShapeDist is then returned.
        Parameters
        ----------
        mapping: LigandAtomMapping
        _rd_shape_dist_func: rdShapeHelpers.ShapeTanimotoDist, optional
            a rdkit function for calculating a shape distance. rdkit so far provides:
            * ShapeProtrudeDist
            * ShapeTanimotoDist
            * ShapeTverskyIndex (however this needs to be curried for alpha and beta parameters)
        _gridSpacing: float, optional
            spacing of the grid for shape calculation, defaults to 0.5
        _vdwScale: float, optional
            range of vdw interactions, defaults to 0.8
        _ignoreHs: bool, optional
            shall Hs be counted as well?, defaults to True
        _maxLayers: int, optional
            number of bits per gridpoints, defaults to -1
        _stepSize: int, optional
            step size between layers, defaults to 0.25

        Returns
        -------
        float
            the ratio of the mol shape distances and the mapped mol shape distances.

        """

        molA = mapping.componentA.to_rdkit()
        molB = mapping.componentB.to_rdkit()
        mapped_atomIDs = mapping.componentA_to_componentB

        #Build mapped mols
        #print("build mols")
        em_A = AllChem.RWMol(molA)
        em_B = AllChem.RWMol(molB)

        i=0
        for atom in molA.GetAtoms():
            if(not atom.GetIdx() in mapped_atomIDs.keys()):
                em_A.RemoveAtom(atom.GetIdx()-i)
                i+=1

        i=0
        for atom in molB.GetAtoms():
            if(not atom.GetIdx() in mapped_atomIDs.values()):
                em_B.RemoveAtom(atom.GetIdx()-i)
                i+=1

        mapping_molA = em_A.GetMol()
        mapping_molB = em_B.GetMol()
        self.mapping_mols = [mapping_molA, mapping_molB]

        #calculate metrics
        #print("calc metrics")
        _shape_dist_func = lambda molA, molB: _rd_shape_dist_func(molA, molB, gridSpacing=_gridSpacing,
                                                                  vdwScale=_vdwScale, ignoreHs=_ignoreHs,
                                                                  maxLayers=_maxLayers, stepSize=_stepSize)
        mol_shape_dist = _shape_dist_func(molA, molB)
        mapped_shape_dist =_shape_dist_func(mapping_molA, mapping_molB)
        shape_dist_ratio = mapped_shape_dist/mol_shape_dist if(mapped_shape_dist != 0 and mol_shape_dist != 0) else 0.0
        print("\tratio", np.round(shape_dist_ratio, 2), "\n\tmol shape D", mol_shape_dist, "\n\tmapped mol shape D", mapped_shape_dist)
        return shape_dist_ratio


class default_kartograf_scorer(_abstract_ligand_atom_mapping_scorer):

    def __init__(self):
        self.scorers = [mapping_ratio_of_mapped_atoms(), mapping_shape_distance()]

    def get_score(self, mapping: LigandAtomMapping) -> float:
        """
            Under Development

        Parameters
        ----------
        mapping : AtomMapping
            ligandAtomMapping to be scored

        Returns
        -------
        float
            normalized score
        """
        d = (scorer(mapping) for scorer in self.scorers)
        score = math.prod(d)
        return score
