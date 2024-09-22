import pytest
import os
import shutil
from pathlib import Path
import mdtraj as md
from openmmdl.openmmdl_simulation.post_md_conversions import MdtrajConverter, MDAnalysisConversion, MDPostProcessingHandler

# Define test data directory and test files
test_data_directory = Path("openmmdl/tests/data/in")
pdb_file = "0_unk_hoh.pdb"
dcd_file = "trajectory.dcd"
ligand_name = "UNK"

@pytest.fixture
def setup_test_environment(tmpdir):
    """Fixture to set up a temporary test environment."""
    original_cwd = os.getcwd()

    # Copy files needed for the test into the temporary directory
    shutil.copy(test_data_directory / pdb_file, original_cwd)
    shutil.copy(test_data_directory / dcd_file, original_cwd)

def test_mdtraj_conversion(setup_test_environment):
    """Test for the MdtrajConverter class."""
    converter = MdtrajConverter(pdb_file, mdtraj_output=["dcd", "xtc", "pdb", "gro"])

    # Run the conversion
    converter.mdtraj_conversion()

    # Check if the output files were created
    assert Path("centered_old_coordinates.dcd").is_file()
    assert Path("centered_old_coordinates.xtc").is_file()
    assert Path("centered_old_coordinates_top.pdb").is_file()
    assert Path("centered_old_coordinates_top.gro").is_file()

def test_mdanalysis_conversion(setup_test_environment):
    """Test for the MDAnalysisConversion class."""
    post_mdtraj_pdb_file = "centered_old_coordinates_top.pdb"
    post_mdtraj_dcd_file = "centered_old_coordinates.dcd"
    converter = MDAnalysisConversion(
        post_mdtraj_pdb_file,
        post_mdtraj_dcd_file,
        mda_output=["pdb", "dcd", "gro", "xtc"],
        output_selection="mda_prot_lig_all",
        ligand_name=ligand_name,
    )

    # Run the conversion
    converter.mdanalysis_conversion()

    # Check if the output files were created
    assert Path("centered_traj.dcd").is_file()
    assert Path("centered_traj_unaligned.dcd").is_file()
    assert Path("centered_top.pdb").is_file()
    assert Path("prot_lig_traj.dcd").is_file()
    assert Path("prot_lig_traj_unaligned.dcd").is_file()
    assert Path("prot_lig_top.pdb").is_file()
    assert Path("centered_traj.xtc").is_file()
    assert Path("centered_traj_unaligned.xtc").is_file()
    assert Path("centered_top.gro").is_file()
    assert Path("prot_lig_traj.xtc").is_file()
    assert Path("prot_lig_traj_unaligned.xtc").is_file()
    assert Path("prot_lig_top.gro").is_file()
