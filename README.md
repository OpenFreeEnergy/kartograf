<p align="center">
    <picture align="center">
      <source media="(prefers-color-scheme: dark)" srcset="https://github.com/OpenFreeEnergy/kartograf/blob/0a2ecb10f1c5ede3842fd9b92ffb232ad726426f/docs/_static/img/kartograf_logo_style5.png">
      <source media="(prefers-color-scheme: light)" srcset="https://github.com/OpenFreeEnergy/kartograf/blob/0a2ecb10f1c5ede3842fd9b92ffb232ad726426f/docs/_static/img/kartograf_logo_style4.png">
      <img alt="Kartografs fancy logo" src="https://github.com/OpenFreeEnergy/kartograf/blob/0a2ecb10f1c5ede3842fd9b92ffb232ad726426f/docs/_static/img/kartograf_logo_style4.png" width=35% >
    </picture>
</p>


Kartograf: A Geometry-Based Atom Mapper
==================================

[//]: # (Badges)
[![Logo](https://img.shields.io/badge/OSMF-OpenFreeEnergy-%23002f4a)](https://openfree.energy/)
[![build](https://github.com/OpenFreeEnergy/kartograf/actions/workflows/ci.yaml/badge.svg)](https://github.com/OpenFreeEnergy/kartograf/actions/workflows/ci.yaml)
[![coverage](https://codecov.io/gh/OpenFreeEnergy/kartograf/branch/main/graph/badge.svg)](https://codecov.io/gh/OpenFreeEnergy/kartograf)
[![Documentation Status](https://readthedocs.org/projects/kartograf/badge/?version=latest)](https://kartograf.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10566716.svg)](https://doi.org/10.5281/zenodo.10566716)
[![Conda Install](https://img.shields.io/badge/Conda%20install---c%20conda--forge%20kartograf-009384)](https://anaconda.org/conda-forge/kartograf)

Kartograf is a package for atom mappings focussing on 3D geometries.
This package can be for example be used to generate hybrid topology systems, where an atom mapping is required to determine the core region of the approach.
But of course there exist also other use cases for this package.
The atom mapper takes two set of coordinates of molecules as input.
Optionally those set of coordinates can be aligned onto each other, checkout the `atom_aligner` module functions 
of Kartograf that offer a shape alignment implementation and a MCS-skeleton alignment.
The `atom_mapper` can be used to generate the 3D geometry focused atom mapping, the algorithm is described in the related publication of Kartograf (see reference).
Additionally, rule based filter functions can be provided to demap atoms, that do not fulfill the desired criteria, see `filters`.
Several mapping scoring metrics are provided, that evaluate geometric properties of your mapping, from `atom_mapping_scorer`, which might be useful for checking quality of your mappings.
Finally, there is a visualization function `display_mappings_3d` that can be used to check out the mappings with a jupyter notebook widget.

Checkout our article on Kartograf in the Journal of Chemical Theory and Computation: [*Kartograf: A Geometrically Accurate Atom Mapper for Hybrid-Topology Relative Free Energy Calculations* - Benjamin Ries*, Irfan Alibay, David W. H. Swenson, Hannah M. Baumann, Michael M. Henry, James R. B. Eastwood, and Richard J. Gowers](https://doi.org/10.1021/acs.jctc.3c01206).

You can find a preprint on [ChemRxiv](https://doi.org/10.26434/chemrxiv-2023-0n1pq).

Try our interactive demo: [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/OpenFreeEnergy/kartograf/blob/main/examples/kartograf_example.ipynb)


## Usage

```python3
from rdkit import Chem
from kartograf.atom_aligner import align_mol_shape
from kartograf.atom_mapping_scorer import MappingRMSDScorer
from kartograf import KartografAtomMapper, SmallMoleculeComponent

# Preprocessing from Smiles - Here you can add your Input!
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

# Score Mapping
rmsd_scorer = MappingRMSDScorer()
score = rmsd_scorer(mapping=kartograf_mapping)
print(f"RMSD Score: {score}")

kartograf_mapping
```
![](docs/_static/img/alignment_benz_ol.png)

## Installation

### Latest release
Kartograf can be installed via the package following package managers:

#### `pip` (PyPI)

```shell
pip install kartograf
```

#### `conda` (conda-forge)

```shell
conda install -c conda-forge kartograf
```

Kartograf can be used via the OpenFE environment like:

```python
from openfe.setup.atom_mapping import kartograf
```

### Development version
The developing setup of Kartograf works like this:

```shell
git clone https://github.com/OpenFreeEnergy/kartograf.git

cd kartograf
mamba env create -f environment.yml

mamba activate kartograf
pip install -e .

```

## License
This library is made available under the MIT open source license.

## Authors

The OpenFE development team.

