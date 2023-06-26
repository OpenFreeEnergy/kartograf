from openfe.utils.visualization_3D import view_mapping_3d
from ipywidgets import interact, widgets


def display_mapping_edges(mappingSet):
    def display_edge(index):
        print("MolA: " + mappingSet[index].componentA.name)
        print("MolB: " + mappingSet[index].componentB.name)
        print("Mapping MolA->MolB:", mappingSet[index].componentA_to_componentB)
        if (hasattr(mappingSet[index], "lscore")):
            print("Lomap Score: ", mappingSet[index].lscore, end="\t")
        if (hasattr(mappingSet[index], "kscore")):
            print("Karto Score: ", mappingSet[index].kscore, end="\n")
        else:
            print()
        view = view_mapping_3d(mappingSet[index], spheres=True, show_atomIDs=True)  # shift=(0.1, 0, 0))
        view.show()

    slider = widgets.IntSlider(tooltip="select mapping",
                               description=str(len(mappingSet)) + " mappings",
                               min=0,
                               max=len(mappingSet) - 1,
                               step=1, value=0)

    nextButton = widgets.Button(
        tooltip='next structure',
        icon='fa-caret-right'
    )

    def increment(fu):
        if (slider.value == slider.max):
            slider.value = 0
        else:
            slider.value += 1

    nextButton.on_click(increment)

    previousButton = widgets.Button(
        tooltip='previous structure',
        icon='fa-caret-left'
    )

    def decrement(fu):
        if (slider.value == 0):
            slider.value = slider.max
        else:
            slider.value -= 1

    previousButton.on_click(decrement)

    hbox = widgets.HBox([previousButton, nextButton, slider])
    inter = widgets.interactive_output(display_edge, {'index': slider})
    vbox = widgets.VBox([hbox, inter])
    return vbox