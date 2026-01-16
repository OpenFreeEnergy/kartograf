# This code is part of kartograf and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/kartograf

from copy import deepcopy
from importlib.resources import files

import pytest
from gufe import SmallMoleculeComponent

from kartograf import KartografAtomMapper
from kartograf.atom_mapper import filter_atoms_h_only_h_mapped
from kartograf.atom_mapper import filter_whole_rings_only
from kartograf.filters.element_change import filter_hybridization_changes
from kartograf.filters.ring_changes import filter_hybridization_rings


def check_mapping_vs_expected(mapping, expected_mapping) -> None:
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
def test_mapping_naphtalene_benzene(naphtalene_benzene_molecules, naphtalene_benzene_mapping) -> None:
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


def test_mapping_naphtalene_benzene_mst(naphtalene_benzene_molecules, naphtalene_benzene_mapping) -> None:
    """
    Test mapping of naphtalene to benzene.
    """
    from kartograf.atom_mapper import mapping_algorithm

    expected_mapping = naphtalene_benzene_mapping.componentA_to_componentB
    geom_mapper = KartografAtomMapper(
        atom_max_distance=0.95,
        atom_map_hydrogens=True,
        map_hydrogens_on_hydrogens_only=False,
        _mapping_algorithm=mapping_algorithm.minimal_spanning_tree,
    )

    geom_mapping = next(
        geom_mapper.suggest_mappings(
            naphtalene_benzene_molecules[0],
            naphtalene_benzene_molecules[1],
        )
    )

    check_mapping_vs_expected(geom_mapping, expected_mapping)


def test_mapping_noMap_algo() -> None:
    """
    Test mapping of naphtalene to benzene.
    """
    with pytest.raises(ValueError) as exc:
        KartografAtomMapper(
            atom_max_distance=0.95,
            atom_map_hydrogens=True,
            map_hydrogens_on_hydrogens_only=False,
            _mapping_algorithm=None,
        )

    assert "Mapping algorithm not implemented or unknown (options: MST or LSA). got key: None" in str(exc.value)


# Check parameters/Filters
def test_mapping_naphtalene_benzene_noHs(naphtalene_benzene_molecules) -> None:
    """
    Test mapping of naphtalene to benzene without H-atoms.
    """
    expected_mapping = {10: 7, 11: 8, 12: 9, 13: 10, 0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5}
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


def test_mapping_naphtalene_benzene_noHs_add_filter(naphtalene_benzene_molecules) -> None:
    """
    Test mapping of naphtalene to benzene without H-atoms added as
     additional filter.
    """
    expected_mapping = {10: 7, 11: 8, 12: 9, 13: 10, 0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5}
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


def test_stereo_mapping(stereco_chem_molecules, stereo_chem_mapping) -> None:
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
    def test_to_from_dict_cycle(self) -> None:
        m = KartografAtomMapper()

        m_dict = m.to_dict()

        m2 = KartografAtomMapper.from_dict(m_dict)

        assert m == m2

    @pytest.mark.parametrize(
        "mhoho,mermo",
        [
            (True, True),
            (True, False),
            (False, True),
            (False, False),
        ],
    )
    def test_check_filters(self, mhoho, mermo) -> None:
        m = KartografAtomMapper(
            map_hydrogens_on_hydrogens_only=mhoho,
            map_exact_ring_matches_only=mermo,
        )

        m2 = KartografAtomMapper.from_dict(m.to_dict())

        assert m._filter_funcs == m2._filter_funcs

    def test_custom_filters(self) -> None:
        def even_prime_filter(mapping):
            mapping.pop(2, None)

            return mapping

        m = KartografAtomMapper(
            additional_mapping_filter_functions=[even_prime_filter],
        )

        m2 = KartografAtomMapper.from_dict(m.to_dict())

        assert len(m._filter_funcs) == len(m2._filter_funcs)

        test_args = {1: 2, 2: 1}
        assert m._filter_funcs[0](test_args) == m2._filter_funcs[0](test_args)


def test_filter_property() -> None:
    mapper = KartografAtomMapper(map_hydrogens_on_hydrogens_only=False)
    first_filters = deepcopy(mapper._filter_funcs)

    mapper.map_hydrogens_on_hydrogens_only = True
    second_filters = deepcopy(mapper._filter_funcs)

    mapper.map_hydrogens_on_hydrogens_only = False
    third_filters = deepcopy(mapper._filter_funcs)

    assert len(first_filters) == len(third_filters)
    assert len(first_filters) == len(second_filters) - 1

    assert filter_atoms_h_only_h_mapped in second_filters
    assert filter_atoms_h_only_h_mapped not in first_filters
    assert filter_atoms_h_only_h_mapped not in first_filters


# Check non params
def test_mapping_rdmols(naphtalene_benzene_molecules, naphtalene_benzene_mapping) -> None:
    """
    Test mapping of naphtalene to benzene.
    """
    naphtalene_benzene_mapping.componentA_to_componentB
    geom_mapper = KartografAtomMapper(
        atom_max_distance=0.95,
        atom_map_hydrogens=True,
        map_hydrogens_on_hydrogens_only=False,
    )

    mols = [naphtalene_benzene_molecules[0], naphtalene_benzene_molecules[1]]

    geom_mapper.suggest_mapping_from_rdmols(
        mols[0].to_rdkit(),
        mols[1].to_rdkit(),
        masked_atoms_molA=None,
        masked_atoms_molB=None,
        pre_mapped_atoms=None,
    )


def test_ring_matches_property() -> None:
    """
    Test ring property changes.
    """
    geom_mapper = KartografAtomMapper(atom_max_distance=0.95, map_exact_ring_matches_only=False)

    geom_mapper.map_exact_ring_matches_only = True
    print([f == filter_whole_rings_only for f in geom_mapper._filter_funcs])
    assert any(f == filter_whole_rings_only for f in geom_mapper._filter_funcs)

    geom_mapper.map_exact_ring_matches_only = False
    assert all(f != filter_whole_rings_only for f in geom_mapper._filter_funcs)

    geom_mapper = KartografAtomMapper(atom_max_distance=0.95, map_exact_ring_matches_only=True)
    assert any(f == filter_whole_rings_only for f in geom_mapper._filter_funcs)


def test_split_multimeric_component() -> None:
    """
    Test splitting chains of 2wtk multimer is done correctly
    """
    from gufe import ProteinComponent
    from openmm.app import PDBFile

    input_pdb = str(files("kartograf.tests.data") / "2wtk_trimer_with_extra_mols.pdb")
    pdb = PDBFile(input_pdb)
    omm_topology = pdb.topology
    # Create the data structure we want to compare to: list[Dict]
    omm_data = []
    # It happens that number of components is the number of chains for this pdb, but doesn't have to
    expected_n_comps = len(list(omm_topology.chains()))
    for chain in omm_topology.chains():
        omm_data.append({"residues": len(list(chain.residues())), "atoms": len(list(chain.atoms()))})

    protein_comp = ProteinComponent.from_pdb_file(input_pdb)
    chain_comps = KartografAtomMapper._split_component_molecules(protein_comp)

    assert len(chain_comps) == expected_n_comps, f"Expected {expected_n_comps} chain components."
    for idx, rdmol in enumerate(chain_comps):
        # TODO: Use MonomerInfo from rdkit to get number of residues. Seems tedious.
        # Make sure the number of atoms in the chain is the same
        n_atoms = rdmol.GetNumAtoms()
        expected_n_atoms = omm_data[idx]["atoms"]
        assert n_atoms == expected_n_atoms, f"Expected {expected_n_atoms}. Received {n_atoms}."


def test_mapping_multimer_components(trimer_2wtk_component, trimer_2wtk_mutated_component) -> None:
    """
    Test we can properly map ProteinComponents generated from 2wtk trimers.

    The final/target component is the same original component but with a ALA-76-TYR mutation.
    """
    from gufe import ProteinComponent

    mapper = KartografAtomMapper(atom_map_hydrogens=True)
    mapping = next(mapper.suggest_mappings(trimer_2wtk_component, trimer_2wtk_mutated_component))
    # It comes from ALA to TYR mutation, n mapped atoms must be 21
    n_atoms_comp_a = trimer_2wtk_component.to_rdkit().GetNumAtoms()
    expected_unique_initial = 1
    expected_mapped_atoms = n_atoms_comp_a - expected_unique_initial
    mapped_atoms = len(mapping.componentA_to_componentB)
    assert mapped_atoms == expected_mapped_atoms, (
        f"Mapped atoms do not match. Expected {expected_mapped_atoms}, received {mapped_atoms}."
    )
    # We expect the unique atoms in initial/ALA to be only 1 hydrogen
    unique_initial = len(list(mapping.componentA_unique))
    assert unique_initial == expected_unique_initial, "Unique atoms in initial molecule do not match."
    # We expect the unique atoms in final/TYR to be 12 atoms
    expected_unique_final = 12
    unique_final = len(list(mapping.componentB_unique))
    assert unique_final == expected_unique_final, "Unique atoms in final molecule do not match."
    # make sure the types and objects have not changed
    assert isinstance(mapping.componentA, ProteinComponent)
    assert isinstance(mapping.componentB, ProteinComponent)
    assert mapping.componentA is trimer_2wtk_component
    assert mapping.componentB is trimer_2wtk_mutated_component


def test_atom_mapping_different_component_types(trimer_2wtk_component, naphtalene_benzene_molecules) -> None:
    """Make sure an error is rasied if we try and create a mapping between two different component types."""
    mapper = KartografAtomMapper()

    with pytest.raises(ValueError, match="were not of the same type, please check the inputs."):
        next(mapper.suggest_mappings(trimer_2wtk_component, naphtalene_benzene_molecules[0]))


def test_atom_mapping_different_number_of_sub_components(trimer_2wtk_component, naphtalene_benzene_molecules) -> None:
    """
    Make sure an error is raised if we get a different number of disconected components in the two molecules
    we want to map.
    """
    mapper = KartografAtomMapper()

    # convert to be the same type to avoid the type check error
    trimer_smc = SmallMoleculeComponent.from_rdkit(trimer_2wtk_component.to_rdkit(), trimer_2wtk_component.name)
    with pytest.raises(
        RuntimeError,
        match="ontain a different number of sub components and so no mapping could be created",
    ):
        next(mapper.suggest_mappings(trimer_smc, naphtalene_benzene_molecules[0]))


@pytest.mark.parametrize(
    "allow_partial_fused_rings, expected_mapping",
    [
        pytest.param(
            True,
            {
                12: 20,
                13: 21,
                14: 22,
                15: 23,
                16: 24,
                17: 19,
                20: 26,
                21: 25,
                0: 6,
                1: 7,
                2: 8,
                3: 9,
                4: 10,
                5: 5,
                6: 4,
                7: 3,
                8: 2,
                9: 13,
                10: 12,
                11: 11,
            },
            id="Allow partial",
        ),
        pytest.param(
            False,
            {12: 20, 13: 21, 14: 22, 15: 23, 16: 24, 0: 6, 1: 7, 2: 8, 3: 9, 4: 10, 5: 5},
            id="Remove partial",
        ),
    ],
)
def test_partial_fused_rings(fused_ring_mols, allow_partial_fused_rings, expected_mapping) -> None:
    """Make sure that partial mappings of fused rings are correctly handled depending on if the flag is set."""
    geom_mapper = KartografAtomMapper(allow_partial_fused_rings=allow_partial_fused_rings)
    geom_mapping = next(
        geom_mapper.suggest_mappings(
            fused_ring_mols[0],
            fused_ring_mols[1],
        )
    )
    check_mapping_vs_expected(mapping=geom_mapping, expected_mapping=expected_mapping)


def test_hybridization_and_ring_breaks(shp2_hybridization_ligands) -> None:
    """
    Make sure rings are not broken when map_exact_ring_matches_only=True and a custom filter is used.
    """
    mapper = KartografAtomMapper(
        map_exact_ring_matches_only=True,
        additional_mapping_filter_functions=[filter_hybridization_changes],
    )
    mapping = next(mapper.suggest_mappings(shp2_hybridization_ligands[0], shp2_hybridization_ligands[1]))
    # check the whole rings are mapped
    filtered_mapping = filter_whole_rings_only(
        mapping.componentA.to_rdkit(),
        mapping.componentB.to_rdkit(),
        mapping.componentA_to_componentB,
    )
    # make sure there was no change in the mapping
    assert filtered_mapping == mapping.componentA_to_componentB


def test_ring_hybridization_with_non_ring_atoms(shp2_hybridization_ligands) -> None:
    """
    Make sure this filter does not fail on non-ring atoms see
    <https://github.com/OpenFreeEnergy/kartograf/issues/62>
    """
    mapper = KartografAtomMapper(additional_mapping_filter_functions=[filter_hybridization_rings])
    mapping = next(mapper.suggest_mappings(shp2_hybridization_ligands[0], shp2_hybridization_ligands[1]))
    # make sure we have some mapping between the atoms
    assert mapping.componentA_to_componentB


@pytest.mark.parametrize(
    "edge, allow_broken, expected_mapping",
    [
        pytest.param(
            ("47", "46"),
            False,
            {
                29: 32,
                30: 31,
                31: 33,
                32: 34,
                33: 38,
                34: 37,
                36: 35,
                37: 36,
                38: 39,
                39: 40,
                40: 41,
                41: 42,
                42: 43,
                43: 44,
                0: 2,
                1: 3,
                2: 0,
                3: 1,
                4: 4,
                5: 5,
                6: 6,
                7: 7,
                8: 8,
                9: 9,
                10: 14,
                11: 10,
                12: 13,
                13: 25,
                14: 11,
                15: 12,
                16: 26,
                17: 15,
                18: 16,
                19: 17,
                20: 18,
                21: 19,
                22: 20,
                23: 21,
                24: 23,
                25: 22,
                26: 24,
            },
            id="47->46",
        ),
        pytest.param(
            ("48", "46"),
            False,
            {
                29: 31,
                30: 32,
                31: 33,
                32: 34,
                33: 37,
                36: 35,
                37: 38,
                38: 39,
                39: 40,
                40: 41,
                41: 42,
                42: 43,
                43: 44,
                0: 0,
                1: 1,
                2: 2,
                3: 3,
                4: 4,
                5: 5,
                6: 6,
                7: 7,
                8: 8,
                9: 9,
                10: 10,
                11: 13,
                12: 25,
                13: 26,
                14: 12,
                15: 11,
                16: 14,
                17: 15,
                18: 16,
                19: 17,
                20: 18,
                21: 19,
                22: 20,
                23: 21,
                24: 23,
                25: 22,
                26: 24,
            },
            id="48->46",
        ),
        pytest.param(
            ("47", "46"),
            True,
            {
                29: 32,
                30: 31,
                31: 33,
                32: 34,
                33: 38,
                34: 37,
                36: 35,
                37: 36,
                38: 39,
                39: 40,
                40: 41,
                41: 42,
                42: 43,
                43: 44,
                0: 2,
                1: 3,
                2: 0,
                3: 1,
                4: 4,
                5: 5,
                6: 6,
                7: 7,
                8: 8,
                9: 9,
                10: 14,
                11: 10,
                12: 13,
                13: 25,
                14: 11,
                15: 12,
                16: 26,
                17: 15,
                18: 16,
                19: 17,
                20: 18,
                21: 19,
                22: 20,
                23: 21,
                24: 23,
                25: 22,
                26: 24,
                27: 28,
                28: 27,
            },
            id="47->46 allow broken",
        ),
    ],
)
def test_bond_break_transforms(pfkfb3_ligands, edge, allow_broken, expected_mapping) -> None:
    """
    Make sure that bond breaking transformations are always filtered, see
    <https://github.com/OpenFreeEnergy/kartograf/issues/88>
    """
    ligand_a = pfkfb3_ligands[edge[0]]
    ligand_b = pfkfb3_ligands[edge[1]]

    mapper = KartografAtomMapper(atom_map_hydrogens=True, allow_bond_breaks=allow_broken)

    mapping = next(mapper.suggest_mappings(ligand_a, ligand_b))
    check_mapping_vs_expected(mapping=mapping, expected_mapping=expected_mapping)
