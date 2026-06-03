====================
Kartograf Change Log
====================

.. current developments

v2.0
====================

**Added:**

* Added documentation for molecule alignment algorithms, including guidance on choosing between skeleton-based MCS alignment and shape-based Open3DAlign alignment (`PR #141 <https://github.com/OpenFreeEnergy/kartograf/issues/141>`_).
* Added a ``kartograf.__version__`` version attribute (PR#132)

**Changed:**

* We now use https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.cdist.html#scipy.spatial.distance.cdist internally to calculate the distance matrix which improves performance.
* **BREAKING CHANGE** The default for ``KartografAtomMapper``'s ``map_hydrogens_on_hydrogens_only`` option is now ``True``.
 Previously it was ``False``.
* **BREAKING CHANGE** We moved the visualization widget to ``gufe``.
  Instead of ``from kartograf.utils.mapping_visualization_widget import display_mappings_3d`` use ``from gufe.visualization.mapping_visualization import display_mappings_3d``

**Deprecated:**

* We no longer distribute the package on PyPI.
  Since all of the dependences of the package do not exist on PyPI, users would have to install the dependences manually.
  We have decided that all new kartograf releases will be solely distributed via conda-forge.
  We have uploaded a placeholder version of the package that when imported, directs users to install the package from conda-forge instead of PyPI.

**Fixed:**

* Clarified atom-mapping warnings when configured filters remove all geometrically found candidate atom pairs (`issue #51 <https://github.com/OpenFreeEnergy/kartograf/issues/51>`_).



v1.2.0
====================

**Added:**

* A flag to retrive the raw kartograf mappings before filtering for bond breaking mappings `PR#90 <https://github.com/OpenFreeEnergy/kartograf/pull/90>`_

**Fixed:**

* Fixed bond breaking transformations from being suggested `PR#89 <https://github.com/OpenFreeEnergy/kartograf/pull/89>`_



v1.1.0
====================

**Added:**

* Added changelog `PR#77 <https://github.com/OpenFreeEnergy/kartograf/pull/77>`_
* Support mapping multi chain components `PR#47 <https://github.com/OpenFreeEnergy/kartograf/pull/47>`_
* Option to not break fused rings when creating mappings `PR#56 <https://github.com/OpenFreeEnergy/kartograf/pull/54>`_
* Added a citation.cff `PR#45 <https://github.com/OpenFreeEnergy/kartograf/pull/45>`_

**Changed:**

* Additional filters are now applied before defaults `PR#64 <https://github.com/OpenFreeEnergy/kartograf/pull/64>`_
* Docs clean up `PR#66 <https://github.com/OpenFreeEnergy/kartograf/pull/66>`_

**Removed:**

* Removed unused scorer `PR#54 <https://github.com/OpenFreeEnergy/kartograf/pull/54>`_

**Fixed:**

* Ring hybridization filter `PR#65 <https://github.com/OpenFreeEnergy/kartograf/pull/65>`_
