# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/kartograf

import math
import numpy as np
from rdkit import Chem
from scipy.spatial import ConvexHull

from gufe import LigandAtomMapping
from gufe.mapping import AtomMapping

import logging
log = logging.getLogger(__name__)

#simple metrics
eukli = lambda x, y: np.sqrt(np.sum(np.square(y - x)))
rms_func = lambda x: np.sqrt(np.mean(np.square(x)))


def mappings_rmsd(mapping: LigandAtomMapping)->float:
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
    for i,j in molA_to_molB.items():
        diff = molPosB[j] - molPosA[i]
        diffs.append(diff)
    
    diffs = np.array(diffs)
    rmsd_map_diff = np.round(np.sqrt(np.sum(diffs**2)),3)
    return rmsd_map_diff


def mapping_mapped_nonMapped_AtomRMS(mapping: AtomMapping):
    raise NotImplementedError()


def mapping_area_ratio(mapping: LigandAtomMapping):
    molA = mapping.componentA.to_rdkit()
    molB = mapping.componentB.to_rdkit()
    molA_to_molB = mapping.componentA_to_componentB
    
    mapping_molA = sorted(list(molA_to_molB.keys()))
    mapping_molB = sorted(list(molA_to_molB.values()))

    if(len(mapping_molA)<4):
        raise ValueError("Mapping is too small to calculate convex hull")

    complete_molA = ConvexHull(molA.GetConformer().GetPositions()).area
    map_molA = ConvexHull(molA.GetConformer().GetPositions()[mapping_molA]).area
    complete_molB = ConvexHull(molB.GetConformer().GetPositions()).area
    map_molB = ConvexHull(molB.GetConformer().GetPositions()[mapping_molB]).area

    ratios = np.array([map_molA/complete_molA, map_molB/complete_molB])
    avg_map = np.mean(ratios)
    
    return 1-avg_map


def mapping_volume_ratio(mapping: LigandAtomMapping)->float:
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
    
    mapping_molA = sorted(list(molA_to_molB.keys()))
    mapping_molB = sorted(list(molA_to_molB.values()))

    if(len(mapping_molA)<4):
        return np.inf
        raise ValueError("Mapping is too small to calculate convex hull")

    complete_molA = ConvexHull(molA.GetConformer().GetPositions()).volume
    map_molA = ConvexHull(molA.GetConformer().GetPositions()[mapping_molA]).volume
    complete_molB = ConvexHull(molB.GetConformer().GetPositions()).volume
    map_molB = ConvexHull(molB.GetConformer().GetPositions()[mapping_molB]).volume

    ratios = np.array([map_molA/complete_molA, map_molB/complete_molB])
    avg_map = np.mean(ratios)
    
    return 1-avg_map


def number_of_mapped_atoms(mapping: LigandAtomMapping)->float:
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
    
    if( len(molB.GetAtoms()) > larger_nAtoms):
        larger_nAtoms = len(molB.GetAtoms())
        
    return len(molA_to_molB)/larger_nAtoms


def default_geometric_score(mapping: LigandAtomMapping)->float:
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
    d = (
        # number_of_mapped_atoms(mapping),
        # mapping_rms(mapping),
        mapping_volume_ratio(mapping),
    )

    score = math.prod(d)
    return score