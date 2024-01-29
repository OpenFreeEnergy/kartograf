import pytest
from rdkit import Chem

from kartograf import filters


@pytest.mark.parametrize('reverse', [False, True])
def test_atoms_H_only_H_mapped(reverse):
    # ethane to propane, hydrogen from ethane mapped to carbon
    m1 = Chem.AddHs(Chem.MolFromSmiles('CC'))
    m2 = Chem.AddHs(Chem.MolFromSmiles('CCC'))
    mapping = {0: 0, 1: 1,
               2: 2,  # this is a H->C, should get removed
               5: 3,  # this is a H->H, should remain
               }
    ref = {0: 0, 1: 1, 5: 3}
    if reverse:
        m1, m2 = m2, m1
        mapping = {v: k for k, v in mapping.items()}
        ref = {v: k for k, v in ref.items()}

    newmapping = filters.filter_atoms_h_only_h_mapped(m1, m2, mapping)

    assert newmapping == ref


@pytest.mark.parametrize('reverse', [False, True])
def test_element_change(reverse):
    # benzene to pyridine, has heteroatom change
    # will result in non-whole ring, but that isn't the job of this filter
    m1 = Chem.MolFromSmiles('c1ccccn1')
    m2 = Chem.MolFromSmiles('c1ccccc1')
    if reverse:
        m1, m2 = m2, m1

    mapping = {i: i for i in range(6)}

    newmapping = filters.filter_element_changes(m1, m2, mapping)

    assert newmapping == {i: i for i in range(5)}

@pytest.mark.parametrize('reverse', [False, True])
def test_element_hybridization_change(reverse):
    # benzene to pyridine, has heteroatom change
    # will result in non-whole ring, but that isn't the job of this filter
    m1 = Chem.MolFromSmiles('CCC')
    m2 = Chem.MolFromSmiles('CC=C')
    if reverse:
        m1, m2 = m2, m1

    mapping = {i: i for i in range(3)}

    newmapping = filters.filter_hybridization_changes(m1, m2, mapping)

    assert newmapping == {i: i for i in range(1)}
