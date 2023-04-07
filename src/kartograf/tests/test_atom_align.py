# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/kartograf

import pytest

from kartograf.atom_align import align_mol_sceletons

from .conf import stereco_chem_molecules


def test_stereo_align(stereco_chem_molecules):
    """
    Currently a smoke test
    """
    molA, molB = stereco_chem_molecules
    aligned_molA = align_mol_sceletons(molA, molB)
