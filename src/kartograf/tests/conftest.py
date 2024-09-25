# This code is part of kartograf and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/kartograf

import pytest
from importlib.resources import files
from rdkit import Chem
from kartograf.atom_aligner import align_mol_skeletons
from gufe import SmallMoleculeComponent, LigandAtomMapping, ProteinComponent


def mol_from_smiles(smiles: str):
    rdmols = [Chem.MolFromSmiles(s) for s in smiles]
    rdmols = [Chem.AddHs(m, addCoords=True) for m in rdmols]
    [
        Chem.rdDistGeom.EmbedMolecule(m, useRandomCoords=False, randomSeed=0)
        for m in rdmols
    ]
    return rdmols


"""
    The Stereochemistry simple Problem
"""


def stereo_chem_mols():
    smiles = [
        "C[C@H](F)Br",
        "C[C@@H](F)Br",
    ]

    rdmols = mol_from_smiles(smiles)
    Chem.rdMolAlign.AlignMol(rdmols[0], rdmols[1])
    mols = [SmallMoleculeComponent(m) for m in rdmols]
    return mols


@pytest.fixture(scope="session")
def stereco_chem_molecules():
    return stereo_chem_mols()


@pytest.fixture(scope="session")
def stereo_chem_mapping():
    mols = stereo_chem_mols()
    expected_mapping = {2: 7, 4: 4, 5: 5, 6: 6, 7: 2, 0: 0, 1: 1, 3: 3}

    return LigandAtomMapping(mols[0], mols[1], expected_mapping)


"""
    napthalene_benzene
"""


def naphtalene_benzene_mols():
    smi_napthalene = "c12ccccc1cccc2"
    smi_benzene = "c1ccccc1"

    rdmols = mol_from_smiles([smi_napthalene, smi_benzene])
    mols = [SmallMoleculeComponent(m) for m in rdmols]
    amol = align_mol_skeletons(mol=mols[1], ref_mol=mols[0])
    mols = [mols[0], amol]

    return mols


@pytest.fixture(scope="session")
def naphtalene_benzene_molecules():
    return naphtalene_benzene_mols()


@pytest.fixture(scope="session")
def naphtalene_benzene_mapping():
    mols = naphtalene_benzene_mols()
    expected_mapping = {
        10: 7,
        11: 8,
        12: 9,
        13: 10,
        0: 0,
        1: 1,
        2: 2,
        3: 3,
        4: 4,
        5: 5,
    }

    return LigandAtomMapping(mols[0], mols[1], expected_mapping)


def benzene_mol():
    smi_benzene = "c1ccccc1"
    rdmols = mol_from_smiles([smi_benzene])
    return SmallMoleculeComponent(rdmols[0])


@pytest.fixture(scope="session")
def benzene_benzene_mapping():
    mol = benzene_mol()
    expected_mapping = {6: 6, 7: 7, 8: 8, 9: 9, 10: 10, 11: 11, 0: 0, 1: 1,
                        2: 2, 3: 3, 4: 4, 5: 5}
    return LigandAtomMapping(mol, mol, expected_mapping)


@pytest.fixture(scope="session")
def benzene_benzene_empty_mapping():
    mol = benzene_mol()
    expected_mapping = {}
    return LigandAtomMapping(mol, mol, expected_mapping)


@pytest.fixture(scope="session")
def trimer_2wtk_component():
    """Protein component for modelled trimer of protein with PDB ID 2wtk"""
    input_pdb = files("kartograf.tests.data").joinpath("2wtk_trimer_with_extra_mols.pdb")
    protein_comp = ProteinComponent.from_pdb_file(str(input_pdb))
    return protein_comp


@pytest.fixture(scope="session")
def trimer_2wtk_mutated_component():
    """Protein component obtained by applying ALA-76-TYR mutation to trimer of 2wtk,
    to residue in chain 'C'."""
    input_pdb = files("kartograf.tests.data").joinpath("2wtk_mutated_with_extra_mols.pdb")
    protein_comp = ProteinComponent.from_pdb_file(str(input_pdb))
    return protein_comp
