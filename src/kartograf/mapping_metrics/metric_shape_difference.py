# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/kartograf

import logging
from typing import Tuple
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdShapeHelpers

from gufe.mapping import AtomMapping

from mapping_metrics._abstract_scorer import _AbstractAtomMappingScorer

log = logging.getLogger(__name__)

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
        elif(mol_shape_dist == 0): # identical mols
            return 0.0
        else:
            shape_dist_ratio = (1*mapped_shape_dist+2*mol_shape_dist)/3
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
