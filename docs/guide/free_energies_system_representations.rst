=======================================
Application in Free Energy Calculations
=======================================

In recent literature many ways showed up on how a system can be represented
during a free energy calculation. Here is a small Taxonomy that tries to
capture them and show benefits and disadvantages.

.. image:: ../_static/img/Topology_types.png

The hybrid topology approach and the single topology approach depend heavily
on atom mappings The mappings are used to find the shared core region of the
molecules, such that the atoms part of this region can be represented as one.
The sigle topology approach tries to maximize here the number of mapped
atoms, where the hybrid topology approach only mapps one shared region of
the molecules and represents the remaining atom coordinates independet of
each other.

The usage of Kartograf's atom mapper for this application can be found in the
turorial!