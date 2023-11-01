import os
import shutil
import pytest
from unittest.mock import mock_open, patch
from openmmdl.openmmdl_simulation.scripts.cleaning_procedures import cleanup, create_directory_if_not_exists, post_md_file_movement, copy_file, create_directory_if_not_exists, organize_files

@pytest.fixture
def test_protein_name():
    return "test_protein"

@pytest.fixture
def test_directory_path():
    return "test_directory"

def test_cleanup(test_protein_name):
    # Create a dummy file to be removed
    with open(f'output_{test_protein_name}', 'w') as dummy_file:
        dummy_file.write("Dummy content")

    # Call the cleanup function
    cleanup(test_protein_name)

    # Check if the file has been removed
    assert not os.path.exists(f'output_{test_protein_name}')

def test_create_directory_if_not_exists(test_directory_path):
    # Create a test directory
    create_directory_if_not_exists(test_directory_path)

    # Check if the directory exists
    assert os.path.exists(test_directory_path)

    # Call the function again, it should not raise an error
    create_directory_if_not_exists(test_directory_path)

    # Cleanup: Remove the test directory
    shutil.rmtree(test_directory_path)
    assert not os.path.exists(test_directory_path)


@patch("os.path.exists")
@patch("shutil.copy")
def test_copy_file(mock_copy, mock_exists):
    
    src = "source_file.txt"
    dest = "destination_directory"

    # Mock the os.path.exists to return True, indicating the source file exists
    mock_exists.return_value = True

    # Call the copy_file function
    copy_file(src, dest)

    # Check that os.path.exists was called with the source file
    mock_exists.assert_called_with(src)

    # Check that shutil.copy was called with the source file and destination directory
    mock_copy.assert_called_with(src, dest)

# Mock the os.path.exists and os.rename functions
@patch("os.path.exists")
@patch("os.rename")
def test_organize_files(mock_rename, mock_exists):
    source = ["file1.txt", "file2.txt", "file3.txt"]
    destination = "destination_directory"

    # Mock os.path.exists to return True for all source files
    mock_exists.side_effect = [True] * len(source)

    # Call the organize_files function
    organize_files(source, destination)

    # Print the calls made to os.rename
    for call in mock_rename.call_args_list:
        print(call)


# Define some sample file paths
protein_name = "sample_protein.pdb"
prmtop = "sample_topology.prmtop"
inpcrd = "sample_coordinates.inpcrd"
ligand = "sample_ligand.pdb"

# Define the fixtures
@pytest.fixture
def mock_organize_files():
    with patch("openmmdl.openmmdl_simulation.scripts.cleaning_procedures.organize_files") as mock:
        yield mock

@pytest.fixture
def mock_copy_file():
    with patch("openmmdl.openmmdl_simulation.scripts.cleaning_procedures.copy_file") as mock:
        yield mock

@pytest.fixture
def mock_create_directory():
    with patch("openmmdl.openmmdl_simulation.scripts.cleaning_procedures.create_directory_if_not_exists") as mock:
        yield mock

def test_post_md_file_movement(
    mock_organize_files, mock_copy_file, mock_create_directory
):
    # Call the function with sample file paths
    post_md_file_movement(protein_name, prmtop, inpcrd, ligand)

    # Check that create_directory_if_not_exists is called with the correct directories
    mock_create_directory.assert_called_with("Input_Files")
    mock_create_directory.assert_called_with("MD_Files/Pre_MD")
    mock_create_directory.assert_called_with("MD_Files/Minimization_Equilibration")
    mock_create_directory.assert_called_with("MD_Files/MD_Output")
    mock_create_directory.assert_called_with("MD_Postprocessing")
    mock_create_directory.assert_called_with("Final_Output")
    mock_create_directory.assert_called_with("Final_Output/All_Atoms")
    mock_create_directory.assert_called_with("Final_Output/Prot_Lig")
    mock_create_directory.assert_called_with("Checkpoints")

    # Check that copy_file is called with the correct files and destinations
    mock_copy_file.assert_called_with(ligand, "Final_Output/All_Atoms")
    mock_copy_file.assert_called_with(ligand, "Final_Output/Prot_Lig")
    mock_copy_file.assert_called_with(protein_name, "Input_Files")
    mock_copy_file.assert_called_with(prmtop, "Input_Files")
    mock_copy_file.assert_called_with(inpcrd, "Input_Files")
    mock_copy_file.assert_called_with(ligand, "Input_Files")

    # Check that organize_files is called with the correct source and destination paths
    expected_source_files = [
        f"{prefix}{protein_name}"
        for prefix in ["prepared_no_solvent_", "solvent_padding_", "solvent_absolute_", "membrane_"]
    ]
    mock_organize_files.assert_called_with(expected_source_files, "MD_Files/Pre_MD")

    expected_source_files = [
        f"{prefix}{protein_name}"
        for prefix in ["Energyminimization_", "Equilibration_"]
    ]
    mock_organize_files.assert_called_with(expected_source_files, "MD_Files/Minimization_Equilibration")

    mock_organize_files.assert_called_with([f"output_{protein_name}", "trajectory.dcd"], "MD_Files/MD_Output")
    mock_organize_files.assert_called_with(
        [
            "centered_old_coordinates_top.pdb",
            "centered_old_coordinates.dcd",
            "centered_old_coordinates_top.gro",
            "centered_old_coordinates.xtc",
        ],
        "MD_Postprocessing",
    )
    mock_organize_files.assert_called_with(
        [
            "centered_top.pdb",
            "centered_traj.dcd",
            "centered_top.gro",
            "centered_traj.xtc",
        ],
        "Final_Output/All_Atoms",
    )
    mock_organize_files.assert_called_with(
        [
            "prot_lig_top.pdb",
            "prot_lig_traj.dcd",
            "prot_lig_top.gro",
            "prot_lig_traj.xtc",
        ],
        "Final_Output/Prot_Lig",
    )
    mock_organize_files.assert_called_with(
        ["checkpoint.chk", "10x_checkpoint.chk", "100x_checkpoint.chk"], "Checkpoints"
    )

# Run the tests
if __name__ == "__main__":
    pytest.main()
