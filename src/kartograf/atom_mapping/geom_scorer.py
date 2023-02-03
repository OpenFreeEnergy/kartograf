# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/kartograf

import math
import numpy as np
from rdkit import Chem

from gufe.mapping import AtomMapping

#simple metrics
eukli = lambda x, y: np.sqrt(np.sum(np.square(y - x)))
rms_func = lambda x: np.sqrt(np.mean(np.square(x)))


def mol_rms(mapping: AtomMapping):
    pass

def mapping_rms(mapping: AtomMapping, _c=1*10**-9):
    molA = mapping.componentA.to_rdkit()
    molB = mapping.componentB.to_rdkit()
    
    posA = molA.GetConformer().GetPositions()
    posB = molB.GetConformer().GetPositions()
    
    atom_map = mapping.componentA_to_componentB
    d = [eukli(np.array(posA[pA]), np.array(posB[pB])) for pA, pB in atom_map.items()]
    max_d = np.max(d)+_c
    rms = rms_func(d)/rms_func(len(d)*max_d)
    
    return rms

def mapping_mapped_nonMapped_AtomRMS(mapping: AtomMapping):
    pass

def mapping_convexVolume_overlap(mapping: AtomMapping):
    pass


def number_of_mapped_atoms(mapping: AtomMapping):
    larger_nAtoms = len(mapping.componentA.to_rdkit().GetAtoms()) 
    if( len(mapping.componentB.to_rdkit().GetAtoms())> larger_nAtoms):
        larger_nAtoms =  len(mapping.componentB.to_rdkit().GetAtoms())
    return len(mapping.componentA_to_componentB)/larger_nAtoms

def geometric_score(mapping: AtomMapping):
    
    d = (
        #number_of_mapped_atoms(mapping),
        mapping_rms(mapping),
    )
    print(d)
    score = math.prod(d)
    return score

    
    