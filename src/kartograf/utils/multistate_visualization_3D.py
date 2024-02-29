'''
3D visualizaiton of multistate mapping.

'''
import numpy as np
from numpy.typing import NDArray
from typing import Union, Optional, Iterable

from rdkit import Chem
from rdkit.Geometry.rdGeometry import Point3D

from matplotlib import pyplot as plt
from matplotlib.colors import rgb2hex

try:
    import py3Dmol
except ImportError:
    pass    # Don't throw  error, will happen later

from gufe import Component
from gufe.mapping import AtomMapping

from kartograf.utils.optional_imports import requires_package



def _translate(mol, shift:Union[tuple[float, float, float], NDArray[
    np.float64]]):
    """
        shifts the molecule by the shift vector

    Parameters
    ----------
    mol : Chem.Mol
        rdkit mol that get shifted
    shift : tuple[float, float, float]
        shift vector

    Returns
    -------
    Chem.Mol
        shifted Molecule (copy of original one)
    """
    mol = Chem.Mol(mol)
    conf = mol.GetConformer()
    for i, atom in enumerate(mol.GetAtoms()):
        x, y, z = conf.GetAtomPosition(i)
        point = Point3D(x + shift[0], y + shift[1], z + shift[2])
        conf.SetAtomPosition(i, point)
    return mol



@requires_package("py3Dmol")
def view_components_3d(mols: Iterable[Component],
                     style: Optional[str] ="stick",
                     shift: Optional[tuple[float, float, float]] = None,
                     view: py3Dmol.view = None
                     ) -> py3Dmol.view:
    """visualize multiple component coordinates in one interactive view.
    It helps to understand how the components are aligned in the system to each other.

    py3Dmol is an optional dependency, it can be installed with:
        pip install py3Dmol

    Parameters
    ----------
    mols : Iterable[ExplicitMoleculeComponent]
        collection of components
    style : Optional[str], optional
        py3Dmol style, by default "stick"
    shift : tuple of floats, optional
        Amount to i*shift each mols_i in order to allow inspection of them in heavy overlap cases.
    view : py3Dmol, optional
        Allows to pass an already existing view, by default None

    Returns
    -------
    py3Dmol.view
        view containing all component coordinates
    """

    if(view is None):
        view = py3Dmol.view(width=600, height=600)

    for i, component in enumerate(mols):
        mol = Chem.Mol(component.to_rdkit())
        if(shift is not None):
            tmp_shift = np.array(shift, dtype=np.float64)*i
            mol = _translate(mol, tmp_shift)

        view.addModel(Chem.MolToMolBlock(mol))

    view.setStyle({style: {}})

    view.zoomTo()
    return view

def visualize_multistate_mappings_3D(components, multi_state_mapping):
    view = view_components_3d(components)
    comp_dict = {c.name: c for c in components}
    cmap = plt.cm.get_cmap("tab20", len(multi_state_mapping))

    for i, am in enumerate(multi_state_mapping):
        color = rgb2hex(cmap(i))

        for lig_name, atom_id in am.items():
            mol = comp_dict[lig_name].to_rdkit()
            p1 = mol.GetConformer().GetAtomPosition(atom_id)

            view.addSphere(
                {
                    "center": {"x": p1.x, "y": p1.y, "z": p1.z},
                    "radius": 0.6,
                    "color": color,
                    "alpha": 0.8,
                }
            )
    return view
