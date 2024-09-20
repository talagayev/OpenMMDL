import os
import shutil
import pytest
from pathlib import Path
from unittest.mock import patch
from openmmdl.openmmdl_simulation.cleaning_procedures import FileManager, Cleanup, PostMDProcessor


@pytest.fixture
def test_protein_name():
    return "test_protein"


@pytest.fixture
def test_directory_path():
    return "test_directory"


def test_cleanup(test_protein_name):
    # Create dummy files to be removed
    with open(f"output_{test_protein_name}", "w") as dummy_file:
        dummy_file.write("Dummy content")
    with open("centered_old_coordinates.pdb", "w") as dummy_file:
        dummy_file.write("Dummy content")
    with open("centered_old_coordinates.dcd", "w") as dummy_file:
        dummy_file.write("Dummy content")

    # Call the cleanup function
    Cleanup.cleanup_files(test_protein_name)

    # Check if the files have been removed
    assert not os.path.exists(f"output_{test_protein_name}")
    assert not os.path.exists("centered_old_coordinates.pdb")
    assert not os.path.exists("centered_old_coordinates.dcd")


def test_create_directory_if_not_exists(test_directory_path):
    # Create a test directory
    FileManager.create_directory_if_not_exists(test_directory_path)

    # Check if the directory exists
    assert os.path.exists(test_directory_path)

    # Call the function again, it should not raise an error
    FileManager.create_directory_if_not_exists(test_directory_path)

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
    FileManager.copy_file(src, dest)

    # Check that os.path.exists was called with the source file
    mock_exists.assert_called_with(src)

    # Check that shutil.copy was called with the source file and destination directory
    mock_copy.assert_called_with(src, dest)


@patch("os.path.exists")
@patch("os.rename")
def test_organize_files(mock_rename, mock_exists):
    source = ["file1.txt", "file2.txt", "file3.txt"]
    destination = "destination_directory"

    # Mock os.path.exists to return True for all source files
    mock_exists.side_effect = [True] * len(source)

    # Call the organize_files function
    FileManager.organize_files(source, destination)

    # Check the os.rename calls
    for call in mock_rename.call_args_list:
        print(call)
    assert mock_rename.call_count == len(source)


def test_post_md_file_movement():
    # Create dummy files and directories for testing
    test_data_directory = Path("openmmdl/tests/data/in")

    ligand = test_data_directory / "CVV.sdf"
    protein_name = test_data_directory / "6b73.pdb"
    prmtop = test_data_directory / "6b73.prmtop"
    inpcrd = test_data_directory / "6b73.inpcrd"
    pre_md_file = test_data_directory / "prepared_no_solvent_6b73.pdb"
    equilibration_file = test_data_directory / "Equilibration_6b73.pdb"
    output_file = test_data_directory / "output_6b73.pdb"
    old_file = test_data_directory / "centered_old_coordinates_top.pdb"
    prot_lig_top = test_data_directory / "prot_lig_top.pdb"
    checkpoint_file = test_data_directory / "checkpoint.chk"

    # Copy files from test_data to the current working directory (like the original pytest might have done)
    cwd = Path.cwd()  # Get current working directory
    shutil.copy(ligand, cwd / "CVV.sdf")
    shutil.copy(protein_name, cwd / "6b73.pdb")
    shutil.copy(prmtop, cwd / "6b73.prmtop")
    shutil.copy(inpcrd, cwd / "6b73.inpcrd")
    shutil.copy(pre_md_file, cwd / "prepared_no_solvent_6b73.pdb")
    shutil.copy(equilibration_file, cwd / "Equilibration_6b73.pdb")
    shutil.copy(old_file, cwd / "centered_old_coordinates_top.pdb")
    shutil.copy(output_file, cwd / "output_6b73.pdb")
    shutil.copy(prot_lig_top, cwd / "prot_lig_top.pdb")
    shutil.copy(checkpoint_file, cwd / "checkpoint.chk")

    # Call the post_md_file_movement function from PostMDProcessor
    PostMDProcessor.post_md_file_movement(
        protein_name="6b73.pdb",
        prmtop="6b73.prmtop",
        inpcrd="6b73.inpcrd",
        ligands="CVV.sdf",
    )

    # Verify the files were moved/copied to the correct locations in MD_Files
    md_files_dir = cwd / "MD_Files"
    md_post_files_dir = cwd / "MD_Postprocessing"
    assert os.path.exists(md_files_dir / "Pre_MD" / "prepared_no_solvent_6b73.pdb")
    assert os.path.exists(md_files_dir / "Minimization_Equilibration" / "Equilibration_6b73.pdb")
    assert os.path.exists(md_post_files_dir / "centered_old_coordinates_top.pdb")
    
    # Check other directories for the copied files
    input_files_dir = cwd / "Input_Files"
    final_output_dir = cwd / "Final_Output"
    checkpoints_dir = cwd / "Checkpoints"
    
    assert os.path.exists(input_files_dir / "6b73.pdb")
    assert os.path.exists(input_files_dir / "6b73.prmtop")
    assert os.path.exists(final_output_dir / "Prot_Lig" / "CVV.sdf")
    assert os.path.exists(final_output_dir / "Prot_Lig" / "prot_lig_top.pdb")
    assert os.path.exists(checkpoints_dir / "checkpoint.chk")

    # Cleanup: Remove created test directories and files
    shutil.rmtree("MD_Files")
    shutil.rmtree("Input_Files")
    shutil.rmtree("Final_Output")
    shutil.rmtree("Checkpoints")

    # Clean up files copied to the current directory
    os.remove("CVV.sdf")
    os.remove("6b73.pdb")
    os.remove("6b73.prmtop")
    os.remove("6b73.inpcrd")


if __name__ == "__main__":
    pytest.main()
