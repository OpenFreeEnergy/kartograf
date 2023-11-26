
Customize Atom Mapping Filter Rules
-----------------------------------
.. _custom-filter-label:

Sometimes the default rule based filters of Kartograf's atom mapper
, that don't map ring size changes, ring breaks or flexibility changes in
rings, might not be enough for you. For example you might want to
additionally avoid element changes in your mapping. Then you can customize
the filter rules and employ your own rules. The general signature of the
Filters is::

    def custom_filter(
        molA: Chem.Mol, molB: Chem.Mol, mapping: dict[int, int]
    ) -> dict[int, int]:
        filtered_mapping = {}

        #do something, but always return a dict with the filtered mapping!

        return filtered_mapping

This signature allows you to build your own filter, with any feature you
like, for example the following definition defines the filter for element
changes, by filtering for the atomic number of the rdkit molecules and
compairing them. The return value of a filter is the filtered dictionary::

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

After defining this filter, you only need to plug it in into the atom mapper
like this::

    from kartograf import KartografAtomMapper
    # Build Kartograf Atom Mapper
    mapper = KartografAtomMapper(additional_mapping_filter_functions=[filter_element_changes])

Now you can start building atom mappings without element changes. Note you
can add as many filters here as you like, they will be exectued in order of
their list appearance. The default ring rules of Kartograf can also be turned
of by setting ``map_exact_ring_matches_only=False``, not recommended though.

