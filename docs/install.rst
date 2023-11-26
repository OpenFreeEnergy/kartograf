=====================
Installation
=====================

User Setup
=============

Kartograf is packaged with OpenFE and can be directly imported from there with:

``from openfe.setup.atom_mapping import kartograf``

Alternatively, you can install Kartograf also as a standalone package via pip
or conda:

``pip install kartograf``

``conda install -c conda-forge kartograf``


Developer Setup
================

For Developers, we recommend to setup a Kartograf environment like in the
following example, allowing modification of the code live in the package::

    git clone https://github.com/OpenFreeEnergy/kartograf.git

    cd kartograf
    mamba env create -f environment.yml

    mamba activate kartograf
    pip install -e .

Happy coding! :)
