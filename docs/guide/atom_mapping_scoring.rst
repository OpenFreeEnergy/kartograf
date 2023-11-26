===============================
Atom Mapping Scoring
===============================
Atom mapping scorers try to evaluate the quality of a given mapping. This
however is a non-trivial task as the quality of a mapping, depends on the use
case of the mapping. A simple case example could be, if a program, for any
reason, reindexes a molecule without changing the conformation, an atom mapping
can be used to get an index matching of before and after the program usage.
Usually, in such a case, you want to map as many atoms as possible with a
distance of 0 angstroms to each other.

A less trivial example is the atom mapping case for molecule transformations
like in hybrid topology free energy calculation approaches. In such a case
the atom quality depends on the likelihood of the transformation to converge
well and to give a reasonable free energy. This in a way would 
require the true and calculated result of the calculation and lead to a henn
egg problem. To tackle this problem usually simplistic approaches
are used to estimate the success likelihood of the transformation. Criteria
for rule-based approaches could be element changes, mapped atom displacement,
flexibility changes or polarity changes, appearing atom number; basically,
anything that introduces a difference from one molecule to the other leading
to an increase of the perturbation. Keep in mind, that this does not only
depend on the ligands themselves but also the environment interacting with
those ligands.

Additionally one can start adding parameter-specific scores, that depend on
force fields or other method-related aspects. However this might
theoretically improve the method outcome, it could minders the
transferability of the scorer from one method to another (overfitting).

In Katograf we added some functionality, that can be used as an aspect of an
atom mapping scorer, like the Mapping RSMD Scorer, checking the displacement
of atoms by the scorers or the Volume Ratio Scorer, checking the volume overlap
of the two molecules.
