# This code is part of kartograf and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/kartograf

from ipywidgets.widgets.widget_box import VBox
from py3Dmol import view

try:
    from kartograf.utils.mapping_visualization_widget import display_mapping_3d, display_mappings_3d

    optional_imports = True
except ImportError:
    optional_imports = False


if optional_imports:

    def test_score_mappings_rmsd(stereo_chem_mapping):
        """
        Currently a smoke test
        """
        v = display_mapping_3d(stereo_chem_mapping)
        assert isinstance(v, view)

    def test_view_mapping(stereo_chem_mapping):
        """
        Currently a smoke test
        """
        view = display_mappings_3d([stereo_chem_mapping, stereo_chem_mapping])
        assert isinstance(view, VBox)
