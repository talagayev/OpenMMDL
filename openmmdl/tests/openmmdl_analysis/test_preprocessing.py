import os
import pytest
from Bio import PDB
import numpy as np
import MDAnalysis as mda
from openmmdl.openmmdl_analysis.preprocessing import process_pdb_file, convert_pdb_to_sdf, renumber_atoms_in_residues, replace_atom_type, process_pdb, move_hydrogens_to_end

pdb_file_path = 'openmmdl/tests/data/in/0_unk_hoh.pdb' 

# Define test data paths
TEST_DATA_DIR = "openmmdl/tests/data/in"
INPUT_PDB_FILENAME = os.path.join(TEST_DATA_DIR, "0_unk_hoh.pdb")
OUTPUT_SDF_FILENAME = os.path.join(TEST_DATA_DIR, "lig.sdf")
OUTPUT_PDB_FILENAME = os.path.join(TEST_DATA_DIR, "0_unk_hoh.pdb")


@pytest.fixture
def sample_pdb_data():
    # Provide sample PDB data for testing
    return "ATOM      1  CA  ASP A   1       0.000   0.000   0.000  1.00  0.00           C"


@pytest.fixture
def temp_pdb_file(tmp_path):
    input_pdb_filename = tmp_path / "test_input.pdb"
    # Copy the content of the provided PDB file to the temporary test file
    with open(pdb_file_path, "r") as src_pdb, open(input_pdb_filename, "w") as dest_pdb:
        dest_pdb.write(src_pdb.read())
    return input_pdb_filename

@pytest.fixture
def sample_pdb_data():
    # Provide sample PDB data for testing
    return "ATOM      1  CA  ASP A   1       0.000   0.000   0.000  1.00  0.00           C"

def test_convert_pdb_to_sdf(tmp_path):
    input_pdb_filename = tmp_path / "input.pdb"
    output_sdf_filename = tmp_path / "output.sdf"
    
    # Create a mock PDB file
    input_pdb_filename.write_text("ATOM      1  CA  ASP A   1       0.000   0.000   0.000  1.00  0.00           C")

    convert_pdb_to_sdf(str(input_pdb_filename), str(output_sdf_filename))
    assert output_sdf_filename.exists()

def test_renumber_atoms_in_residues(sample_pdb_data, tmp_path):
    input_pdb_filename = tmp_path / "input.pdb"
    output_pdb_filename = tmp_path / "output.pdb"

    # Create a mock PDB file
    input_pdb_filename.write_text(sample_pdb_data)

    renumber_atoms_in_residues(str(input_pdb_filename), str(output_pdb_filename), 'ASP')
    assert output_pdb_filename.exists()

def test_replace_atom_type(tmp_path):
    input_file = tmp_path / "input.pdb"
    output_file = tmp_path / "output.pdb"

    # Create a mock PDB file
    input_file.write_text("ATOM      1  CA  LIG X   1       0.000   0.000   0.000  1.00  0.00           C")

    # Read input data from the mock input file
    with open(input_file, 'r') as f:
        input_data = f.read()

    # Apply the replace_atom_type function to the input data
    modified_data = replace_atom_type(input_data)

    # Write the modified data to the mock output file
    with open(output_file, 'w') as f:
        f.write(modified_data)

    # Check if the modification was successful
    assert 'LIG C' in output_file.read_text()
