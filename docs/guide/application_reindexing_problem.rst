=======================================
Application for Reindexing Problem
=======================================

Atom indices are important properties of Molecules in any insilico approach.
However it still sometimes happens, that there are pieces of code, that
reindex the atoms of a molecule like in the following example:

.. image:: ../_static/img/reindexing_problem_vis.png

If the molecule has a given conformation, that stays unchanged, Kartograf's
atom mapper can be used to address this problem very easily and precisely!
The result here would be the mapping from input to output molecule: `{1:7, 2:4, 3:3, 4:2, 5:1, 6:6, 7:5}`

You can use the code from the mapping tutorial in order to solve this problem:
:doc:`/tutorial/mapping_tutorial`.
