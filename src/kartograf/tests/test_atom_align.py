# This code is part of kartograf and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/kartograf


from kartograf.atom_aligner import align_mol_shape, align_mol_skeletons


def test_stereo_align_mcs(stereco_chem_molecules) -> None:
    """
    Currently a smoke test
    """
    molA, molB = stereco_chem_molecules
    align_mol_skeletons(molA, molB)


def test_stereo_align_shape(stereco_chem_molecules) -> None:
    """
    Currently a smoke test
    """
    molA, molB = stereco_chem_molecules
    align_mol_shape(molA, molB)
