# This code is part of kartograf and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/kartograf

import pytest

from kartograf.atom_aligner import align_mol_shape, align_mol_skeletons

from .conftest import stereco_chem_molecules


def test_stereo_align_mcs(stereco_chem_molecules):
    """
    Currently a smoke test
    """
    molA, molB = stereco_chem_molecules
    aligned_molA = align_mol_skeletons(molA, molB)


def test_stereo_align_shape(stereco_chem_molecules):
    """
    Currently a smoke test
    """
    molA, molB = stereco_chem_molecules
    aligned_molA = align_mol_shape(molA, molB)
