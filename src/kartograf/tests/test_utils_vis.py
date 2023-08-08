import pytest

from py3Dmol import view
from ipywidgets.widgets.widget_box import VBox

try:
    from kartograf.utils.mapping_visualization_widget import display_mappings_3d, display_mapping_3d
    optional_imports = True
except ImportError:
    optional_imports = False

from .conf import stereo_chem_mapping


if(optional_imports):
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
        view = display_mappings_3d([stereo_chem_mapping,stereo_chem_mapping])
        assert isinstance(view, VBox)
