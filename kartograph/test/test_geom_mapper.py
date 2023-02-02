import pytest

from rdkit import Chem


from kartograph import atom_mapping




@pytest.fixture(scope="session")
def benzols():
    
    m = Chem.MolFromPDBBlock(thrombin_pdb, proximityBonding=False, removeHs=False)

    return pdbinf.assign_pdb_bonds(m, templates=[pdbinf.STANDARD_AA_DOC])


class test_geomety_atom_mapper():
    
    def test_mapping(t)