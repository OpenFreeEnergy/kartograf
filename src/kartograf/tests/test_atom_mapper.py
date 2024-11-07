# This code is part of kartograf and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/kartograf

import pytest
from importlib.resources import files
from kartograf import KartografAtomMapper
from kartograf.atom_mapper import (filter_atoms_h_only_h_mapped,
                                   filter_whole_rings_only)

from .conftest import (
    naphtalene_benzene_molecules,
    naphtalene_benzene_mapping,
    stereco_chem_molecules,
    stereo_chem_mapping,
)
from copy import deepcopy


def check_mapping_vs_expected(mapping, expected_mapping):
    assert len(expected_mapping) == len(mapping.componentA_to_componentB)

    diff = []
    for exp_k, exp_v in expected_mapping.items():
        if exp_k not in mapping.componentA_to_componentB:
            diff.append((exp_k, exp_v, None, "missing mapping for atom"))
        elif exp_v != mapping.componentA_to_componentB[exp_k]:
            diff.append(
                (
                    exp_k,
                    exp_v,
                    mapping.componentA_to_componentB[exp_k],
                    "-",
                    "differently mapped atoms",
                )
            )

    if len(diff) > 0:
        print("Differences: expected key, expected value, actual value")
        print("\n".join(map(lambda x: "\t".join(map(str, x)), diff)))
        raise ValueError("mapping did not match expected mapping!")


# Mapping Algorithm tests
def test_mapping_naphtalene_benzene(
    naphtalene_benzene_molecules, naphtalene_benzene_mapping
):
    """
    Test mapping of naphtalene to benzene.
    """
    expected_mapping = naphtalene_benzene_mapping.componentA_to_componentB
    geom_mapper = KartografAtomMapper(
        atom_max_distance=0.95,
        atom_map_hydrogens=True,
        map_hydrogens_on_hydrogens_only=False,
    )

    geom_mapping = next(
        geom_mapper.suggest_mappings(
            naphtalene_benzene_molecules[0],
            naphtalene_benzene_molecules[1],
        )
    )

    check_mapping_vs_expected(geom_mapping, expected_mapping)


def test_mapping_naphtalene_benzene_mst(
    naphtalene_benzene_molecules, naphtalene_benzene_mapping
):
    """
    Test mapping of naphtalene to benzene.
    """
    from kartograf.atom_mapper import mapping_algorithm
    expected_mapping = naphtalene_benzene_mapping.componentA_to_componentB
    geom_mapper = KartografAtomMapper(
        atom_max_distance=0.95,
        atom_map_hydrogens=True,
        map_hydrogens_on_hydrogens_only=False,
        _mapping_algorithm=mapping_algorithm.minimal_spanning_tree
    )

    geom_mapping = next(
        geom_mapper.suggest_mappings(
            naphtalene_benzene_molecules[0],
            naphtalene_benzene_molecules[1],
        )
    )

    check_mapping_vs_expected(geom_mapping, expected_mapping)


def test_mapping_noMap_algo():
    """
    Test mapping of naphtalene to benzene.
    """
    with pytest.raises(ValueError) as exc:
        geom_mapper = KartografAtomMapper(
            atom_max_distance=0.95,
            atom_map_hydrogens=True,
            map_hydrogens_on_hydrogens_only=False,
            _mapping_algorithm=None
        )

    assert "Mapping algorithm not implemented or unknown (options: MST or " \
           "LSA). got key: None" in str(exc.value)


# Check parameters/Filters
def test_mapping_naphtalene_benzene_noHs(naphtalene_benzene_molecules):
    """
    Test mapping of naphtalene to benzene without H-atoms.
    """
    expected_mapping = {10: 7, 11: 8, 12: 9, 13: 10, 0: 0, 1: 1, 2: 2, 3: 3,
                        4: 4, 5: 5}
    geom_mapper = KartografAtomMapper(
        atom_max_distance=0.95,
        atom_map_hydrogens=True,
        map_hydrogens_on_hydrogens_only=True,
    )

    geom_mapping = next(
        geom_mapper.suggest_mappings(
            naphtalene_benzene_molecules[0],
            naphtalene_benzene_molecules[1],
        )
    )

    check_mapping_vs_expected(geom_mapping, expected_mapping)


def test_mapping_naphtalene_benzene_noHs_add_filter(naphtalene_benzene_molecules):
    """
    Test mapping of naphtalene to benzene without H-atoms added as
     additional filter.
    """
    expected_mapping = {10: 7, 11: 8, 12: 9, 13: 10, 0: 0, 1: 1, 2: 2,
                        3: 3, 4: 4, 5: 5}
    geom_mapper = KartografAtomMapper(
        atom_max_distance=0.95,
        atom_map_hydrogens=True,
        map_hydrogens_on_hydrogens_only=False,
        additional_mapping_filter_functions=[filter_atoms_h_only_h_mapped],
    )

    geom_mapping = next(
        geom_mapper.suggest_mappings(
            naphtalene_benzene_molecules[0],
            naphtalene_benzene_molecules[1],
        )
    )

    check_mapping_vs_expected(geom_mapping, expected_mapping)


def test_stereo_mapping(stereco_chem_molecules, stereo_chem_mapping):
    """
    Test weird Stereochemistry mapping.
    """
    expected_mapping = stereo_chem_mapping.componentA_to_componentB
    geom_mapper = KartografAtomMapper(
        atom_max_distance=0.95, atom_map_hydrogens=True
    )  # mapping_algorithm.minimal_spanning_tree
    geom_mapping = next(
        geom_mapper.suggest_mappings(
            stereco_chem_molecules[0],
            stereco_chem_molecules[1],
        )
    )

    check_mapping_vs_expected(geom_mapping, expected_mapping)


class TestSerialisation:
    def test_to_from_dict_cycle(self):
        m = KartografAtomMapper()

        m_dict = m.to_dict()

        m2 = KartografAtomMapper.from_dict(m_dict)

        assert m == m2

    @pytest.mark.parametrize('mhoho,mermo', [
        (True, True),
        (True, False),
        (False, True),
        (False, False),
    ])
    def test_check_filters(self, mhoho, mermo):
        m = KartografAtomMapper(
            map_hydrogens_on_hydrogens_only=mhoho,
            map_exact_ring_matches_only=mermo,
        )

        m2 = KartografAtomMapper.from_dict(m.to_dict())

        assert m._filter_funcs == m2._filter_funcs

    def test_custom_filters(self):
        def nop_filter(a, b, c):
            return c

        m = KartografAtomMapper(
            additional_mapping_filter_functions=[nop_filter],
        )

        m2 = KartografAtomMapper.from_dict(m.to_dict())

        assert m._filter_funcs == m2._filter_funcs


def test_filter_property():
    mapper = KartografAtomMapper(map_hydrogens_on_hydrogens_only=False)
    first_filters = deepcopy(mapper._filter_funcs)

    mapper.map_hydrogens_on_hydrogens_only = True
    second_filters = deepcopy(mapper._filter_funcs)

    mapper.map_hydrogens_on_hydrogens_only = False
    third_filters = deepcopy(mapper._filter_funcs)

    assert len(first_filters) == len(third_filters)
    assert len(first_filters) == len(second_filters)-1

    assert filter_atoms_h_only_h_mapped in second_filters
    assert filter_atoms_h_only_h_mapped not in first_filters
    assert filter_atoms_h_only_h_mapped not in first_filters


# Check non params
def test_mapping_rdmols(
    naphtalene_benzene_molecules, naphtalene_benzene_mapping
):
    """
    Test mapping of naphtalene to benzene.
    """
    expected_mapping = naphtalene_benzene_mapping.componentA_to_componentB
    geom_mapper = KartografAtomMapper(
        atom_max_distance=0.95,
        atom_map_hydrogens=True,
        map_hydrogens_on_hydrogens_only=False,
    )

    mols = [naphtalene_benzene_molecules[0],
            naphtalene_benzene_molecules[1]]

    m = geom_mapper.suggest_mapping_from_rdmols(mols[0].to_rdkit(),
                                                mols[1].to_rdkit(),
                                                masked_atoms_molA=None,
                                                masked_atoms_molB=None,
                                                pre_mapped_atoms=None)


def test_ring_matches_property():
    """
    Test ring property changes.
    """
    geom_mapper = KartografAtomMapper(
        atom_max_distance=0.95,
        map_exact_ring_matches_only=False
    )

    geom_mapper.map_exact_ring_matches_only = True
    print([f == filter_whole_rings_only for f in geom_mapper._filter_funcs])
    assert any(f == filter_whole_rings_only for f in geom_mapper._filter_funcs)

    geom_mapper.map_exact_ring_matches_only = False
    assert all(f != filter_whole_rings_only for f in geom_mapper._filter_funcs)

    geom_mapper = KartografAtomMapper(
        atom_max_distance=0.95,
        map_exact_ring_matches_only=True
    )
    assert any(f == filter_whole_rings_only for f in geom_mapper._filter_funcs)


def test_split_multimeric_component():
    """
    Test splitting chains of 2wtk multimer is done correctly
    """
    from openmm.app import PDBFile
    from kartograf.atom_mapper import KartografAtomMapper
    from gufe import ProteinComponent

    input_pdb = str(files("kartograf.tests.data").joinpath("2wtk_trimer_with_extra_mols.pdb"))
    pdb = PDBFile(input_pdb)
    omm_topology = pdb.topology
    # Create the data structure we want to compare to: list[Dict]
    omm_data = []
    # It happens that number of components is the number of chains for this pdb, but doesn't have to
    expected_n_comps = len(list(omm_topology.chains()))
    for chain in omm_topology.chains():
        omm_data.append({"residues": len(list(chain.residues())),
                         "atoms": len(list(chain.atoms()))})

    protein_comp = ProteinComponent.from_pdb_file(input_pdb)
    chain_comps = KartografAtomMapper._split_protein_components_molecules(protein_comp)

    assert len(chain_comps) == expected_n_comps, f"Expected {expected_n_comps} chain components."
    for idx, comp in enumerate(chain_comps):
        rdmol = comp.to_rdkit()
        # TODO: Use MonomerInfo from rdkit to get number of residues. Seems tedious.
        # Make sure the number of atoms in the chain is the same
        n_atoms = rdmol.GetNumAtoms()
        expected_n_atoms = omm_data[idx]["atoms"]
        assert n_atoms == expected_n_atoms, f"Expected {expected_n_atoms}. Received {n_atoms}."
    

def test_mapping_multimer_components(trimer_2wtk_component,
                                     trimer_2wtk_mutated_component):
    """
    Test we can properly map ProteinComponents generated from 2wtk trimers.

    The final/target component is the same original component but with a ALA-76-TYR mutation.
    """
    from kartograf import KartografAtomMapper
    mapper = KartografAtomMapper(atom_map_hydrogens=True)
    mapping = next(mapper.suggest_mappings(trimer_2wtk_component,
                                           trimer_2wtk_mutated_component))
    # It comes from ALA to TYR mutation, n mapped atoms must be 21
    n_atoms_comp_a = trimer_2wtk_component.to_rdkit().GetNumAtoms()
    expected_unique_initial = 1
    expected_mapped_atoms = n_atoms_comp_a - expected_unique_initial
    mapped_atoms = len(mapping.componentA_to_componentB)
    assert mapped_atoms == expected_mapped_atoms, \
        f"Mapped atoms do not match. Expected {expected_mapped_atoms}, received {mapped_atoms}."
    # We expect the unique atoms in initial/ALA to be only 1 hydrogen
    unique_initial = len(list(mapping.componentA_unique))
    assert unique_initial == expected_unique_initial, \
        f"Unique atoms in initial molecule do not match."
    # We expect the unique atoms in final/TYR to be 12 atoms
    expected_unique_final = 12
    unique_final = len(list(mapping.componentB_unique))
    assert unique_final == expected_unique_final, f"Unique atoms in final molecule do not match."
