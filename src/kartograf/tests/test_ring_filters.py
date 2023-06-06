import pytest
from rdkit import Chem

from kartograf import filters


@pytest.mark.parametrize('reverse', [False, True])
def test_ringsize_filter(reverse):
    # naphthalene to indole, 6,6 to 6,5
    # should trim out the 5->6 ring
    m1 = Chem.MolFromSmiles('c12ccccc1cccc2')
    m2 = Chem.MolFromSmiles('c12ccccc1C=CN2')
    if reverse:
        m1, m2 = m2, m1

    mapping = {i: i for i in range(9)}

    newmapping = filters.filter_ringsize_changes(m1, m2, mapping)

    assert newmapping == {i: i for i in range(6)}


def test_ringsize_safe():
    m1 = Chem.MolFromSmiles('c1ccccc1')
    m2 = Chem.MolFromSmiles('c1ccccn1')

    mapping = {i: i for i in range(6)}
    newmapping = filters.filter_ringsize_changes(m1, m2, mapping)

    assert newmapping == mapping


@pytest.mark.parametrize('reverse', [False, True])
def test_ringbreaks(reverse):
    # naphthalene to toluene
    # should remove methyl
    m1 = Chem.MolFromSmiles('c12ccccc1cccc2')
    m2 = Chem.MolFromSmiles('c1ccccc1C')

    mapping = {i: i for i in range(7)}

    if reverse:
        m1, m2 = m2, m1

    newmapping = filters.filter_ringbreak_changes(m1, m2, mapping)

    assert newmapping == {i: i for i in range(6)}


def test_ringbreaks_safe():
    m1 = Chem.MolFromSmiles('c1ccccc1')
    m2 = Chem.MolFromSmiles('c1ccccn1')

    mapping = {i: i for i in range(6)}
    newmapping = filters.filter_ringbreak_changes(m1, m2, mapping)

    assert newmapping == mapping


@pytest.mark.parametrize('reverse', [False, True])
def test_whole_rings_only(reverse):
    # benzene to pyridine
    # only first 5 atoms mapped, so none should remain mapped
    m1 = Chem.MolFromSmiles('c1ccccc1')
    m2 = Chem.MolFromSmiles('c1ccccn1')
    if reverse:
        m1, m2 = m2, m1

    mapping = {i: i for i in range(5)}

    newmapping = filters.filter_whole_rings_only(m1, m2, mapping)

    assert newmapping == {}


def test_whole_rings_safe():
    m1 = Chem.MolFromSmiles('c1ccccc1')
    m2 = Chem.MolFromSmiles('c1ccccn1')

    mapping = {i: i for i in range(6)}
    newmapping = filters.filter_whole_rings_only(m1, m2, mapping)

    assert newmapping == mapping
