.. _custom-filter-label:

Customize Atom Mapping Filter Rules
-----------------------------------

Sometimes the default rule-based filters of Kartograf's atom mapper,
which don't map ``ring size changes``, ``ring breaks``, or ``flexibility changes in rings``, might not be enough for you.
For example, you might want to additionally avoid element changes in your mapping. Then you can customize
the filter rules and employ your own rules. The general signature of the filters uses two molecules ``molA`` and
``molB``, that were provided to generate a ``mapping`` linking individual atoms from ``molA`` to ``molB``.
In the function, the ``mapping`` should then be filtered by the implemented rule resulting in a returned
``filtered_mapping``:

.. code-block::

    def custom_filter(
        molA: Chem.Mol, molB: Chem.Mol, mapping: dict[int, int]
    ) -> dict[int, int]:
        filtered_mapping = {}

        #do something, but always return a dict with the filtered mapping!

        return filtered_mapping

This signature allows you to build your own filter, with any feature you
like, for example, the following definition defines the filter for the element
changes, by filtering for the atomic number of the RDKit molecules and
comparing them. The return value of a filter is the filtered dictionary:

.. code-block::

    def filter_element_changes(
        molA: Chem.Mol, molB: Chem.Mol, mapping: dict[int, int]
    ) -> dict[int, int]:
        """Forces a mapping to exclude any alchemical element changes in the core"""
        filtered_mapping = {}

        for i, j in mapping.items():
            if (
                molA.GetAtomWithIdx(i).GetAtomicNum()
                != molB.GetAtomWithIdx(j).GetAtomicNum()
            ):
                continue
            filtered_mapping[i] = j

        return filtered_mapping

After defining this filter, you only need to plug it into the atom mapper
like this:

.. code-block::

    from kartograf import KartografAtomMapper
    # Build Kartograf Atom Mapper
    mapper = KartografAtomMapper(additional_mapping_filter_functions=[filter_element_changes])

Now you can start building atom mappings without element changes. Note you
can add as many filters here as you like, they will be executed in order of
their list appearance. The default ring rules of Kartograf can also be turned
off by setting ``map_exact_ring_matches_only=False``, not recommended though.
