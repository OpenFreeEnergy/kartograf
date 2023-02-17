from rdkit.Chem import rdFMCS
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

from gufe import SmallMoleculeComponent


def align_molB_to_molA_sceletons(molA:SmallMoleculeComponent, 
                                 molB:SmallMoleculeComponent,
                                 verbose:bool = False,
                        ) -> SmallMoleculeComponent:
    mol1b = molA._rdkit
    mol2b = molB._rdkit
        
    p = rdFMCS.MCSParameters()
    p.AtomTyper = rdFMCS.AtomCompare.CompareAny


    res = rdFMCS.FindMCS([mol1b, mol2b],p)

    # convert match to mapping'
    q = Chem.MolFromSmarts(res.smartsString)
    if(verbose): print(q)

    m1_idx = mol1b.GetSubstructMatch(q)
    m2_idx = mol2b.GetSubstructMatch(q)
    if(verbose): print(m1_idx, m2_idx)

    idx_mappings = list(zip(m2_idx,m1_idx))

    rms = AllChem.AlignMol(
        prbMol=mol2b, refMol=mol1b, atomMap=idx_mappings,
    )
    if(verbose): print(rms)

    molB._rdkit = mol2b
    return molB


