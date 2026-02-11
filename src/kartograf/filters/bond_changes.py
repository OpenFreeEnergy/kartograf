# This code is part of kartograf and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/kartograf
import logging

import networkx as nx
from rdkit import Chem

logger = logging.getLogger(__name__)


def filter_bond_breaks(mol_a: Chem.Mol, mol_b: Chem.Mol, mapping: dict[int, int]):
    """
    Remove cases where a bond would be broken during the transformation.
    These are detected by checking that all bonds of ligand A which are in the mapping also form a bond in ligand B,
    if they do not the atoms are removed from the mapping.

    See <https://github.com/OpenFreeEnergy/kartograf/issues/88> for more details.

    Parameters
    ----------
    mol_a
        Molecule at state A with atom indices in the keys of the mapping dict
    mol_b
        Molecule at state B with atom indices in the values of the mapping dict
    mapping
        The mapping between molecule A and B

    Returns
    -------
        A filtered mapping removing any entries which map a broken bond.
    """

    # generate a list of bonds in molecule A and B
    mol_a_bonds = [sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())) for bond in mol_a.GetBonds()]
    mol_b_bonds = [sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())) for bond in mol_b.GetBonds()]
    broken_bonds = []

    # create a graph to get the largest connected mapping if we have a bond break
    graph = nx.Graph()

    # loop over the bonds in molecule A and check they are connected in molecule B via the mapping
    for atom_a, atom_b in mol_a_bonds:
        if atom_a in mapping and atom_b in mapping:
            # the bond is in the mapped section
            # make sure the mapped atoms are also bonded
            mapped_bond = sorted((mapping[atom_a], mapping[atom_b]))
            if mapped_bond in mol_b_bonds:
                # if the bond is conserved add it to the graph
                graph.add_edge(atom_a, atom_b)

            else:
                logger.debug(f"Bond {mapped_bond} is not preserved between molecule A and B possibly a broken bond.")
                broken_bonds.append((atom_a, atom_b))
                continue

    # if there is a broken bond try and find the largest submapping
    if broken_bonds:
        sub_graphs = [subgraph for subgraph in nx.connected_components(graph)]
        # sort by length to get the largest one
        sub_graphs.sort(key=lambda x: len(x), reverse=True)
        largest_graph = sub_graphs[0]
        mapping = {k: v for k, v in mapping.items() if k in largest_graph}

    return mapping
