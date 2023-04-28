<p align="center">
    <img src=".img/Kartograf_logo_boxed_dark_transp.png" class="gh-dark-mode-only"/>
    <img src=".img/Kartograf_logo_boxed_light_transp.png" class="gh-light-mode-only"/>
</p>

Kartograf: A 3D atom graph mapper
==================================

[//]: # (Badges)
[![Logo](https://img.shields.io/badge/OSMF-OpenFreeEnergy-%23002f4a)](https://openfree.energy/)
[![build](https://github.com/OpenFreeEnergy/kartograf/actions/workflows/ci.yaml/badge.svg)](https://github.com/OpenFreeEnergy/kartograf/actions/workflows/ci.yaml)
[![coverage](https://codecov.io/gh/OpenFreeEnergy/kartograf/branch/main/graph/badge.svg)](https://codecov.io/gh/OpenFreeEnergy/kartograf)
[![Documentation Status](https://readthedocs.org/projects/kartograf/badge/?version=latest)](https://kartograf.readthedocs.io/en/latest/?badge=latest)

[![Pip Install](https://img.shields.io/badge/pip%20install-kartograf-d9c4b1)](https://img.shields.io/badge/pip%20install-kartograf-d9c4b1)
[![Conda Install](https://img.shields.io/badge/Conda%20install---c%20conda--forge%20kartograf-009384)](https://img.shields.io/badge/Conda%20install---c%20conda--forge%20kartograf-009384)


Kartograf offers a geometric atom mapper approach, that allows to map a given set of ligand coordinates. (can be used for hybrid topology  RBFE calculations)
This package can be used standalone, or from the OpenFE environment.

**More will be here soon!**

## Usage
```python3
from rdkit import Chem
from kartograf.atom_align import align_mol_shape
from kartograf import KartografAtomMapper, SmallMoleculeComponent

#Preprocessing from Smiles - Here you can add your Input!
# Generate Data: START
smiles = ["c1ccccc1", "c1ccccc1(CO)"]
rdmols = [Chem.MolFromSmiles(s) for s in smiles]
rdmols = [Chem.AddHs(m, addCoords=True) for m in rdmols]
[Chem.rdDistGeom.EmbedMolecule(m, useRandomCoords=False, randomSeed = 0) for m in rdmols]
# Generate Data: END

# Build Small Molecule Components
molA, molB = [SmallMoleculeComponent.from_rdkit(m) for m in rdmols]

# Align the mols first - this might not needed, depends on input.
a_molB = align_mol_shape(molB, ref_mol=molA)


# Build Kartograf Atom Mapper
mapper = KartografAtomMapper(atom_map_hydrogens=True)

# Get Mapping
kartograf_mapping = next(mapper.suggest_mappings(molA, a_molB))

kartograf_mapping
```
![](.img/alignment_benz_ol.png)

## Installation
you can install Kartograf via the package manager of your choice:

```shell
pip install kartograf
```

```shell
conda install -c conda-forge kartograf
```

Or use Kartograf from the OpenFE Environment (soon).

## References



