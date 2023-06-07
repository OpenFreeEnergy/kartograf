# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/kartograf

import numpy as np
from typing import Tuple

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdShapeHelpers

from scipy.spatial import ConvexHull
from scipy import constants as const

from gufe import AtomMapping, SmallMoleculeComponent
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
class _AbstractAtomMappingScorer(abc.ABC):

    def __init__(self):
        pass

    def __call__(self, mapping:AtomMapping, *args, **kwargs)->float:
        return self.get_score(mapping)

    @abc.abstractmethod
    def get_score(self,mapping:AtomMapping, *args, **kwargs)->float:
        """
            the scoring function returns a value between 0 and 1.
            a value close to Zero indicates a small distance, a score close to one indicates a large cost/error.

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

class MappingRMSDScorer(_AbstractAtomMappingScorer):
    def get_rmsd(self, mapping: AtomMapping) -> float:
        """this function calculates the rmsd between the mapped atoms of the two molecules

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
        rmsd_map_diff = np.round(np.sqrt(np.sum(diffs**2)), 3)

        return float(np.round(rmsd_map_diff,5))

    def get_rmsd_p(self, mapping:AtomMapping, accepted_distance_rmsd:float = 0.5, k_hook=1, T=298) -> float:
        """estimate likelihood of this shift by calculating the probability of the rmsd of the mapping with a harmonic oscillator

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
            1 - likelihood of the shift (1-p)
        """
        rmsd = self.get_rmsd(mapping)
        beta = 1000 / (const.k * const.Avogadro * T)
        V = k_hook * (rmsd-accepted_distance_rmsd)
        p = np.exp(-beta * V) if(np.exp(-beta * V) < 1) else 1
        return float(np.round(p,5))

    def get_score(self, mapping, accepted_distance_rmsd:float = 0, k_hook=10**0, T=298)->float:
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
        return float(np.round(1-self.get_rmsd_p(mapping, accepted_distance_rmsd = accepted_distance_rmsd, k_hook= k_hook, T= T),2))


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
        return 1 if(r<0) else np.round(1-r, 2)

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

        if len(mapping_molA) < 4:
            return 0.0
            #raise ValueError("Mapping is too small to calculate convex hull")

        complete_molA = ConvexHull(molA.GetConformer().GetPositions()).volume
        map_molA = ConvexHull(molA.GetConformer().GetPositions()[mapping_molA]).volume
        complete_molB = ConvexHull(molB.GetConformer().GetPositions()).volume
        map_molB = ConvexHull(molB.GetConformer().GetPositions()[mapping_molB]).volume

        ratios = np.array([map_molA / complete_molA, map_molB / complete_molB])
        avg_map_volume_ratio = np.mean(ratios)

        #print("ratios",avg_map_volume_ratio, ratios)
        #print("volumes", map_molA, complete_molA, map_molB, complete_molB)
        #print('ind', mapping_molA, mapping_molB)
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

        return np.round(1-(len(molA_to_molB) / larger_nAtoms), 2)


class _MappingShapeDistanceScorer(_AbstractAtomMappingScorer):
    mapping_mols: Tuple[Chem.Mol, Chem.Mol]

    def __init__(self, _rd_shape_dist_func=rdShapeHelpers.ShapeTanimotoDist, _grid_spacing:float=0.5, _vdw_scale:float=0.8, _ignore_hs:bool=False,
                 _max_layers:int=-1, _step_size:float=0.25):
        """
        This function is using the implemented shape distances in rdkit and applies them to the atom mapping problem.


        Parameters
        ----------
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


        """
        super().__init__()
        self._shape_dist_func = lambda molA, molB: _rd_shape_dist_func(molA, molB, gridSpacing=_grid_spacing,
                                                                       vdwScale=_vdw_scale, ignoreHs=_ignore_hs,
                                                                       maxLayers=_max_layers, stepSize=_step_size)

    def get_score(self, mapping:AtomMapping) ->float:
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
        s = self.get_mapping_shape_distance(mapping)
        return np.round(s,2) if(s<1) else 1.0

    def get_mapped_mols(self, mapping:AtomMapping)->Tuple[Chem.Mol, Chem.Mol  ]:
        """
        reduce the two molecules to an rdmol, representing the mapping.

        Parameters
        ----------
        mapping: AtomMapping
            AtomMapping to be scored

        Returns
        -------
        Chem.Mol, Chem.Mol
            returns the two mapping based molecules starting once from molA and once from molB
        """
        molA = mapping.componentA.to_rdkit()
        molB = mapping.componentB.to_rdkit()
        mapped_atomIDs = mapping.componentA_to_componentB

        if(len(mapped_atomIDs) == 0):
            return [None, None]

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
        return [mapping_molA, mapping_molB]


    def get_rdmol_shape_distance(self, molA:Chem.Mol, molB: Chem.Mol)->float:
        """
        Calculates the shape distance of two rdmols.

        Parameters
        ----------
        molA: Chem.Mol
        molB: Chem.Mol

        Returns
        -------
        float
            raw result of the applied rdkit shape distance.
        """
        if(any([x is None for x in [molA, molB]])):
            if(all([x is None for x in [molA, molB]])):
                return np.inf
            else:
                raise ValueError("One rdkit mol is None! (molA, molB)", [molA, molB])
        return self._shape_dist_func(molA=molA, molB=molB)

    def get_mol_shape_distance(self, molA:SmallMoleculeComponent, molB:SmallMoleculeComponent)->float:
        """
        Calculates the shape distance of two gufe Components.

        Parameters
        ----------
        molA: Chem.Mol
        molB: Chem.Mol

        Returns
        -------
        float
            raw result of the applied rdkit shape distance.
        """
        return self.get_rdmol_shape_distance(molA=molA.to_rdkit(), molB=molB.to_rdkit())

    def get_mapping_mol_shape_distance(self, mapping:AtomMapping)->float:
        """
        Calculates the shape distance of a gufe AtomMapping .

        Parameters
        ----------
        mapping: AtomMapping

        Returns
        -------
        float
            raw result of the applied rdkit shape distance.
        """
        return self.get_mol_shape_distance(molA=mapping.componentA, molB=mapping.componentB)

    def get_mapped_structure_shape_distance(self, mapping:AtomMapping)->float:
        """
        Calculates the shape distance of the mapped parts of the molecules only.

        Parameters
        ----------
        mapping: AtomMapping

        Returns
        -------
        float
            raw result of the applied rdkit shape distance.
        """
        self.mapping_mols = self.get_mapped_mols(mapping=mapping)
        mapped_molA, mapped_molB = self.mapping_mols
        return self.get_rdmol_shape_distance(molA=mapped_molA, molB=mapped_molB)

    def get_mapping_shape_distance(self, mapping:AtomMapping)->float:
        """
            the function builds two mapped mols, each for one molecule. Basically reduces the atoms to the mapped region.
            next it calculates the shape distance of both full molecules and of both mapped mols.
            The ratio mappedMolShapeDist/MolShapeDist is then returned.
        Parameters
        ----------
        mapping: AtomMapping


        Returns
        -------
        float
            the ratio of the mol shape distances and the mapped mol shape distances.

        """

        #calculate metrics
        mol_shape_dist = self.get_mapping_mol_shape_distance(mapping=mapping)
        mapped_shape_dist = self.get_mapped_structure_shape_distance(mapping=mapping)

        if(mapped_shape_dist == np.inf):    # no mapping
            return 1.0
        elif(mol_shape_dist == 0 and mapped_shape_dist == 0): # identical mols
            return 0.0
        else:
            shape_dist_ratio = mapped_shape_dist/mol_shape_dist
            return shape_dist_ratio


class MappingShapeOverlapScorer(_MappingShapeDistanceScorer):
    def __init__(self, _grid_spacing=0.5, _vdw_scale=0.8, _ignore_hs=False, _max_layers=-1, _step_size=0.25):
        """
        This class uses the _MappingShapeDistanceScorer with the settings such, that the overlap of the two molecules are taken into consideration.

        Parameters
        ----------
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
        """

        super().__init__(_rd_shape_dist_func=rdShapeHelpers.ShapeTanimotoDist, _grid_spacing=_grid_spacing,
                         _vdw_scale=_vdw_scale, _ignore_hs=_ignore_hs, _max_layers=_max_layers, _step_size=_step_size)

class MappingShapeMismatchScorer(_MappingShapeDistanceScorer):

    def __init__(self, _grid_spacing=0.5, _vdw_scale=0.8, _ignore_hs=False, _max_layers=-1, _step_size=0.25):
        """
        This class uses the _MappingShapeDistanceScorer with the settings such, that the volume mismatches of the two molecules are taken into consideration.

        Parameters
        ----------
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
        """

        super().__init__(_rd_shape_dist_func=rdShapeHelpers.ShapeProtrudeDist, _grid_spacing=_grid_spacing,
                         _vdw_scale=_vdw_scale, _ignore_hs=_ignore_hs, _max_layers=_max_layers, _step_size=_step_size)

class DefaultKartografScorer(_AbstractAtomMappingScorer):
    """
    Warning this is highly experimental!
    """
    def __init__(self):
        """
        this could be a future scorer using the here derined metrics?

        don't use ;)
        """
        self.scorers = [MappingVolumeRatioScorer(), #MappingRatioMappedAtomsScorer
                        MappingShapeOverlapScorer(),
                        MappingShapeMismatchScorer(),
                        MappingRMSDScorer()]

        self.weights = np.array([1,3,3,3])
        self.weigths = self.weights/sum(self.weights)

    def get_score(self, mapping: AtomMapping) -> float:
        """
            Under Development

        Parameters
        ----------
        mapping : AtomMapping
            AtomMapping to be scored

        Returns
        -------
        float
            normalized score
        """
        log.info("Kartograf Score:")
        scores = []
        for weight, scorer in zip(self.weights, self.scorers):
            s = scorer.get_score(mapping)
            log.info("\t"+scorer.__class__.__name__+"\t"+str(s)+"\tweight: "+str(weight))
            scores.append(s*weight)

        score = np.round(np.mean(scores),2)
        score = score if(score<1.0) else 1.0
        log.info("Result: "+str(score))
        return score
