Comparing Two Mappings
-----------------------
If you want to compare two mappings, maybe from different approaches like
Kartograf and Lomap, you can of course always compare the number of mapped atoms::

    nA = len(mappingA.componentA_to_componentB)
    nB = len(mappingB.componentA_to_componentB)
    nA>nB

This will give you insights how many atoms get mapped by your mapping
approaches. But this will not tell you how similar actually the mapped atoms
are. For that we added the Jaccard Score (jcs) to the Kartograf repository.
The jcs is calculated as follows

.. math:: s_{jaccard} = |A \cap B| / |A \cup B|
   :label: jaccard_score

The score can be used to calculate the atom pair diversity of the mapping A and mapping B::

    # an atom mapping - atom_mapping
    from kartograf.metrics_mapping_comparisons import jaccard_score

    jcs = jaccard_score(mappingA, mappingB)

This will return you a number between `0` and `1`. If your jsc is `1`, it means
your mappings are identical. But if it is `0` you have entirely different
mappings, mapping different atoms from mol A to mol B.
