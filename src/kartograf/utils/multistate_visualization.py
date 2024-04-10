

'''
2D
'''

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from IPython.display import Image, display


def visualize_multistate_mappings_2D(components, multi_state_mapping, ncols=5):
    nrows = len(components) // ncols
    nrows = nrows if (len(components) % ncols == 0) else nrows + 1
    grid_x = ncols
    grid_y = nrows
    d2d = Draw.rdMolDraw2D.MolDraw2DCairo(grid_x * 500, grid_y * 500, 500, 500)

    # squash to 2D
    copies = [Chem.Mol(mol.to_rdkit()) for mol in components]
    for mol in copies:
        AllChem.Compute2DCoords(mol)

    # mol alignments if atom_mapping present
    ref_mol = copies[0]
    for mobile_mol in copies[1:]:
        atomMap = []
        for ms_map in multi_state_mapping:
            atomMap.append((ms_map[mobile_mol.GetProp("_Name")],
                            ms_map[ref_mol.GetProp("_Name")]))

        AllChem.AlignMol(mobile_mol, ref_mol, atomMap=atomMap)

    atom_lists = []
    for c in components:
        lig_maps = []
        for m in multi_state_mapping:
            lig_maps.append(m[c.name])
        atom_lists.append(lig_maps)

    RED = (220 / 255, 50 / 255, 32 / 255, 1.0)
    # standard settings for our visualization
    d2d.drawOptions().useBWAtomPalette()
    d2d.drawOptions().continousHighlight = False
    d2d.drawOptions().setHighlightColour(RED)
    d2d.drawOptions().addAtomIndices = True
    d2d.DrawMolecules(
        copies,
        highlightAtoms=atom_lists,
        # highlightBonds=bonds_list,
        # highlightAtomColors=atom_colors,
        # highlightBondColors=bond_colors,
    )
    d2d.FinishDrawing()

    return Image(d2d.GetDrawingText())



