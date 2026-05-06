import sys
from unittest.mock import MagicMock, patch

import pytest
from plip.basic import config

from openmmdl.openmmdl_analysis.analysis.interactions import InteractionAnalyzer


def test_process_trajectory_prolif_peptide_empty_selection_raises_clean_error(monkeypatch):
    """Wrong peptide chain IDs should raise a clean ValueError, not UnboundLocalError."""
    monkeypatch.setitem(sys.modules, "prolif", MagicMock())

    analyzer = InteractionAnalyzer.__new__(InteractionAnalyzer)
    analyzer.special = None
    analyzer.peptide = "B"
    analyzer.pdb_md = MagicMock()

    empty_ligand_ag = MagicMock()
    empty_ligand_ag.__len__.return_value = 0
    receptor_ag = MagicMock()
    receptor_ag.__len__.return_value = 10
    analyzer.pdb_md.select_atoms.side_effect = [empty_ligand_ag, receptor_ag]

    with pytest.raises(ValueError, match="ProLIF: ligand selection returned 0 atoms: chainID B"):
        analyzer._process_trajectory_prolif()


@patch("openmmdl.openmmdl_analysis.analysis.interactions.PDBComplex")
def test_retrieve_plip_interactions_peptide_sets_config_peptides(mock_pdb_complex, monkeypatch):
    """Peptide PLIP workers should not rely on inherited parent-process config."""
    monkeypatch.setattr(config, "PEPTIDES", [], raising=False)

    analyzer = InteractionAnalyzer.__new__(InteractionAnalyzer)
    analyzer.peptide = "B"

    protlig = mock_pdb_complex.return_value
    protlig.ligands = []

    with pytest.raises(ValueError, match="PLIP did not detect peptide chain 'B' as a ligand"):
        analyzer._retrieve_plip_interactions_peptide("processing_frame_1.pdb")

    assert config.PEPTIDES == ["B"]
    protlig.load_pdb.assert_called_once_with("processing_frame_1.pdb")
