# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/kartograf

from typing import List, Union
from ipywidgets import widgets

from gufe import AtomMapping

try:
    from openfe.utils.visualization_3D import (
        view_mapping_3d as display_mapping_3d,
    )
except ImportError:
    pass  # Don't throw  error, will happen later

from .optional_imports import requires_package


@requires_package("py3Dmol")
def display_mappings_3d(mappingSet: Union[AtomMapping, List[AtomMapping]]) -> widgets.VBox:
    """
    This function is visualizing the provided list of mappings. It shows in the middle an overlay of the coordinates
    of the molecues, and left and right the mapping of the atoms (color of the spheres indicates partners).
    This function is tested in jupyter notebooks.

    Parameters
    ----------
    mappingSet:Union[AtomMapping, List[AtomMapping]]
        a list of atom mappings (gufe.AtomMapping objects)

    Returns
    -------
    widgets.Vbox:
        returns a widget, with the visualization and control elements.

    """
    #Input Parse
    if(isinstance(mappingSet, AtomMapping)):
        mappingSet = [mappingSet]

    #helper for drawing edges
    def display_edge(index):
        print("MolA: " + mappingSet[index].componentA.name)
        print("MolB: " + mappingSet[index].componentB.name)
        print(
            "Mapping MolA->MolB:", mappingSet[index].componentA_to_componentB
        )
        if hasattr(mappingSet[index], "score"):
            print("Mapping Score: ", getattr(mappingSet[index], "score"))
        else:
            print()
        view = display_mapping_3d(
            mappingSet[index], spheres=True, show_atomIDs=True
        )  # shift=(0.1, 0, 0))
        view.show()

    # Int slider for selecting mapping from set
    slider = widgets.IntSlider(
        tooltip="select mapping",
        description=str(len(mappingSet)) + " mappings",
        min=0,
        max=len(mappingSet) - 1,
        step=1,
        value=0,
    )

    # jump one mapping forward/backwards
    nextButton = widgets.Button(
        tooltip="next structure", icon="fa-caret-right"
    )

    def increment(fu):
        if slider.value == slider.max:
            slider.value = 0
        else:
            slider.value += 1

    nextButton.on_click(increment)

    previousButton = widgets.Button(
        tooltip="previous structure", icon="fa-caret-left"
    )

    def decrement(fu):
        if slider.value == 0:
            slider.value = slider.max
        else:
            slider.value -= 1

    previousButton.on_click(decrement)

    #Aligning control elements and visualization
    hbox = widgets.HBox([previousButton, nextButton, slider])
    inter = widgets.interactive_output(display_edge, {"index": slider})
    vbox = widgets.VBox([hbox, inter])

    return vbox
