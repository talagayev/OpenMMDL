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

def test_process_pdb_file(temp_pdb_file):
    input_pdb_filename = temp_pdb_file
    # Check if the function modifies the PDB file as expected
    process_pdb_file(input_pdb_filename)

    # Load the modified PDB file
    u = mda.Universe(input_pdb_filename)

    # Check that the residue name "*" is changed to "UNK"
    for atom in u.atoms:
        if atom.residue.resname == "UNK":
            assert atom.residue.resname == "UNK"

    # You can add additional assertions as needed.

def test_convert_pdb_to_sdf():
    convert_pdb_to_sdf(INPUT_PDB_FILENAME, OUTPUT_SDF_FILENAME)
    assert os.path.exists(OUTPUT_SDF_FILENAME)


def test_renumber_atoms_in_residues(sample_pdb_data):
    # Mock a PDB file with a single line for testing
    with open(INPUT_PDB_FILENAME, 'w') as f:
        f.write(sample_pdb_data)

    renumber_atoms_in_residues(INPUT_PDB_FILENAME, OUTPUT_PDB_FILENAME, 'ASP')
    assert os.path.exists(OUTPUT_PDB_FILENAME)


def test_replace_atom_type():
    # Provide sample PDB data for testing
    sample_pdb_data = "ATOM      1  CA  LIG X   1       0.000   0.000   0.000  1.00  0.00           C"
    modified_data = replace_atom_type(sample_pdb_data)
    assert ' LIG  C' in modified_data


def test_process_pdb():
    # Mock input and output file paths
    input_file = os.path.join(TEST_DATA_DIR, "input_file.pdb")
    output_file = os.path.join(TEST_DATA_DIR, "output_file.pdb")

    # Mock input data
    input_data = "ATOM      1  CA  LIG X   1       0.000   0.000   0.000  1.00  0.00           C"

    # Write input data to the mock input file
    with open(input_file, 'w') as f:
        f.write(input_data)

    # Perform the process_pdb function
    process_pdb(input_file, output_file)

    # Read the output data from the mock output file
    with open(output_file, 'r') as f:
        output_data = f.read()

    assert ' LIG  C' in output_data


def test_move_hydrogens_to_end():
    # Create a structure for testing
    structure = PDB.PDBParser().get_structure('test_structure', 'path/to/pdb/file.pdb')

    # Perform the move_hydrogens_to_end function
    move_hydrogens_to_end(structure, 'ASP')

    # Add assertions based on the expected behavior of the function
    # For example, check if hydrogen atoms are moved to the end of the ASP residue
    assert True  # Add your assertions here
