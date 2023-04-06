import pytest
from rdkit import Chem
from gufe import SmallMoleculeComponent
from kartograf.atom_align import align_mol_sceletons
from gufe import SmallMoleculeComponent, LigandAtomMapping

def mol_from_smiles(smiles: str):
    rdmols = [Chem.MolFromSmiles(s) for s in smiles]
    rdmols = [Chem.AddHs(m, addCoords=True) for m in rdmols]
    [Chem.rdDistGeom.EmbedMolecule(m, useRandomCoords=False, randomSeed = 0) for m in rdmols]
    return rdmols


@pytest.fixture(scope="session")
def stereco_chem_molecules():
    smiles = [
        "C[C@H](F)Br",
        "C[C@@H](F)Br",
    ]

    rdmols = mol_from_smiles(smiles)
    Chem.rdMolAlign.AlignMol(rdmols[0], rdmols[1])
    mols = [SmallMoleculeComponent(m) for m in rdmols]
    return mols

@pytest.fixture(scope="session")
def stereo_mapping():
    smiles = [
        "C[C@H](F)Br",
        "C[C@@H](F)Br",
    ]

    mols = [Chem.MolFromSmiles(s) for s in smiles]
    mols = [Chem.AddHs(m, addCoords=True) for m in mols]
    [Chem.rdDistGeom.EmbedMultipleConfs(m, 1) for m in mols]
    Chem.rdMolAlign.AlignMol(mols[0], mols[1])

    return LigandAtomMapping(
        SmallMoleculeComponent(mols[0]),
        SmallMoleculeComponent(mols[1]),
        {2: 7, 4: 5, 5: 4, 6: 6, 7: 2, 0: 0, 1: 1, 3: 3},
    )

@pytest.fixture(scope="session")
def naphtalene_benzene():
    smi_napthalene = 'c12ccccc1cccc2'
    smi_benzene = 'c1ccccc1'

    rdmols = mol_from_smiles([smi_napthalene, smi_benzene])
    mols = [SmallMoleculeComponent(m) for m in rdmols]
    amol = align_mol_sceletons(mol=mols[1], ref_mol=mols[0])
    mols = [mols[0], amol]

    return mols