import os
import shutil
import pytest
from openmmdl.openmmdl_simulation.scripts.cleaning_procedures import cleanup, create_directory_if_not_exists

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
