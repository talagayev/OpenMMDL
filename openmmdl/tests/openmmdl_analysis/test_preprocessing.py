import os
import pytest
import tempfile
import shutil
from Bio import PDB
import numpy as np
import mdtraj as md
from pathlib import Path
import MDAnalysis as mda
from openmmdl.openmmdl_analysis.preprocessing import *

pdb_file_path = "openmmdl/tests/data/in/0_unk_hoh.pdb"

# Define test data paths
test_data_directory = Path("openmmdl/tests/data/in")
pdb_file = test_data_directory / "0_unk_hoh.pdb"
topology_metal = f"{test_data_directory}/metal_top.pdb"
ligand_resname = "UNK"


# Import the renumber_protein_residues function from your module
from your_module_name import renumber_protein_residues

# Define paths to test PDB files
input_pdb_path = "test_input.pdb"
reference_pdb_path = "test_reference.pdb"
output_pdb_path = "test_output.pdb"  # New output PDB path

# Fixture to create sample input PDB file
@pytest.fixture
def create_input_pdb():
    # Save the sample input PDB content to a file
    input_pdb_content = """
    ATOM      1  N   ALA A   1      11.112  31.224  51.546  1.00  0.00
    ATOM      2  CA  ALA A   1      11.207  31.756  50.180  1.00  0.00
    ATOM      3  C   ALA A   1      12.635  31.978  49.720  1.00  0.00
    ATOM      4  O   ALA A   1      13.349  32.835  50.348  1.00  0.00
    TER
    """
    with open(input_pdb_path, "w") as f:
        f.write(input_pdb_content)
    return input_pdb_path

# Fixture to create sample reference PDB file
@pytest.fixture
def create_reference_pdb():
    # Save the sample reference PDB content to a file with the changed residue number
    reference_pdb_content = """
    ATOM      1  N   ALA A   5      10.112  30.224  50.546  1.00  0.00
    ATOM      2  CA  ALA A   5      10.207  30.756  49.180  1.00  0.00
    ATOM      3  C   ALA A   5      11.635  30.978  48.720  1.00  0.00
    ATOM      4  O   ALA A   5      12.349  31.835  49.348  1.00  0.00
    TER
    """
    with open(reference_pdb_path, "w") as f:
        f.write(reference_pdb_content)
    return reference_pdb_path

# Fixture to create the output PDB path
@pytest.fixture
def create_output_pdb_path():
    return output_pdb_path

# Fixture to clean up test files after testing
@pytest.fixture
def clean_up_files():
    # Code to clean up test files after testing
    yield
    if os.path.exists(input_pdb_path):
        os.remove(input_pdb_path)
    if os.path.exists(reference_pdb_path):
        os.remove(reference_pdb_path)
    if os.path.exists(output_pdb_path):
        os.remove(output_pdb_path)

# Pytest function to test renumber_protein_residues
def test_renumber_protein_residues(create_input_pdb, create_reference_pdb, create_output_pdb_path, clean_up_files):
    # Call the renumber_protein_residues function with the test files
    renumber_protein_residues(create_input_pdb, create_reference_pdb, create_output_pdb_path)

    # Load the output PDB file and check if the renumbering is correct
    output_traj = md.load(output_pdb_path)
    input_traj = md.load(create_input_pdb)
    reference_traj = md.load(create_reference_pdb)

    # Additional assertions based on your specific requirements
    assert output_traj.n_frames == input_traj.n_frames
    assert output_traj.n_atoms == input_traj.n_atoms
    # Add more assertions if needed

    # Check if the residue indices are correctly renumbered
    for frame in range(output_traj.n_frames):
        for chain_id in output_traj.topology.chain_ids:
            output_residues = list(output_traj.topology.residues_in_chain(chain_id))
            input_residues = list(input_traj.topology.residues_in_chain(chain_id))
            reference_residues = list(reference_traj.topology.residues_in_chain(chain_id))

            for i, output_residue in enumerate(output_residues):
                output_residue_index = output_residue.index
                input_residue_index = input_residues[i].index
                reference_residue_index = reference_residues[i].index

                # Add assertions to check if the renumbering is correct
                assert output_residue_index == reference_residue_index - input_residue_index


@pytest.fixture
def sample_pdb_data():
    # Provide sample PDB data for testing
    return """ATOM      1  N   UNK A 454      43.493  48.319  35.835  1.00  0.00      A    N  
ATOM      2  N1  UNK A 454      44.740  47.862  35.697  1.00  0.00      A    N  
ATOM      3  C14 UNK A 454      44.608  46.866  34.829  1.00  0.00      A    C  
ATOM      4  N2  UNK A 454      43.265  46.644  34.450  1.00  0.00      A    N  
ATOM      5  C7  UNK A 454      42.607  47.556  35.077  1.00  0.00      A    C  
ATOM      6  H5  UNK A 454      41.542  47.701  34.954  1.00  0.00      A    H  
ATOM      7  H10 UNK A 454      45.308  46.132  34.453  1.00  0.00      A    H  
ATOM      8  C   UNK A 454      43.168  49.513  36.656  1.00  0.00      A    C  
ATOM      9  C2  UNK A 454      42.743  50.705  35.818  1.00  0.00      A    C  
ATOM     10  C4  UNK A 454      43.545  51.052  34.671  1.00  0.00      A    C  
ATOM     11  C9  UNK A 454      43.171  52.151  33.897  1.00  0.00      A    C  
ATOM     12  C13 UNK A 454      42.090  52.924  34.222  1.00  0.00      A    C  
ATOM     13  C11 UNK A 454      41.393  52.671  35.378  1.00  0.00      A    C  
ATOM     14  C6  UNK A 454      41.793  51.635  36.268  1.00  0.00      A    C  
ATOM     15  H4  UNK A 454      41.220  51.358  37.148  1.00  0.00      A    H  
ATOM     16  H9  UNK A 454      40.518  53.291  35.552  1.00  0.00      A    H  
ATOM     17  C16 UNK A 454      41.790  54.079  33.432  1.00  0.00      A    C  
ATOM     18  N4  UNK A 454      41.594  54.934  32.652  1.00  0.00      A    N  
ATOM     19  H7  UNK A 454      43.694  52.248  32.951  1.00  0.00      A    H  
ATOM     20  H2  UNK A 454      44.333  50.369  34.369  1.00  0.00      A    H  
ATOM     21  H   UNK A 454      44.108  49.790  37.148  1.00  0.00      A    H  
ATOM     22  C1  UNK A 454      42.146  49.054  37.737  1.00  0.00      A    C  
ATOM     23  C5  UNK A 454      42.675  48.761  39.003  1.00  0.00      A    C  
ATOM     24  C10 UNK A 454      41.859  48.278  39.998  1.00  0.00      A    C
ATOM     25  H8  UNK A 454      42.284  48.099  40.981  1.00  0.00      A    H  
ATOM     26  H3  UNK A 454      43.752  48.806  39.135  1.00  0.00      A    H  
ATOM     27  C3  UNK A 454      40.774  48.885  37.463  1.00  0.00      A    C  
ATOM     28  H1  UNK A 454      40.310  49.079  36.500  1.00  0.00      A    H  
ATOM     29  C8  UNK A 454      39.907  48.435  38.509  1.00  0.00      A    C  
ATOM     30  H6  UNK A 454      38.833  48.310  38.406  1.00  0.00      A    H  
ATOM     31  C12 UNK A 454      40.466  48.125  39.823  1.00  0.00      A    C  
ATOM     32  C15 UNK A 454      39.627  47.605  40.833  1.00  0.00      A    C  
ATOM     33  N3  UNK A 454      38.981  47.235  41.740  1.00  0.00      A    N """


@pytest.fixture
def input_pdb_filename(tmp_path):
    input_pdb_filename = tmp_path / "input.pdb"
    with open(input_pdb_filename, "w") as f:
        f.write(
            """ATOM      1  N   SPC A 101      43.493  48.319  35.835  1.00  0.00      A    N  
ATOM      2  O   TIP3 A 102      44.740  47.862  35.697  1.00  0.00      A    O  
ATOM      3  C   *   A 103      44.608  46.866  34.829  1.00  0.00      A    C  
ATOM      4  H   *   A 104      43.265  46.644  34.450  1.00  0.00      A    H  
ATOM      5  O   WAT A 105      42.607  47.556  35.077  1.00  0.00      A    O  
ATOM      6  H1  SPC A 106      41.542  47.701  34.954  1.00  0.00      A    H  
ATOM      7  H2  *   A 107      45.308  46.132  34.453  1.00  0.00      A    H  
ATOM      8  C   T3P A 108      43.168  49.513  36.656  1.00  0.00      A    C  
ATOM      9  O   T4P A 109      42.743  50.705  35.818  1.00  0.00      A    O  
ATOM     10  H   T5P A 110      43.545  51.052  34.671  1.00  0.00      A    H  
ATOM     11  N   *   A 111      43.171  52.151  33.897  1.00  0.00      A    N  
ATOM     12  C   SPC A 112      42.090  52.924  34.222  1.00  0.00      A    C  
ATOM     13  O   *   A 113      41.393  52.671  35.378  1.00  0.00      A    O  
ATOM     14  C   TIP4 A 114      41.793  51.635  36.268  1.00  0.00      A    C  
ATOM     15  O   *   A 115      41.220  51.358  37.148  1.00  0.00      A    O  
ATOM     16  H   *   A 116      40.518  53.291  35.552  1.00  0.00      A    H  
ATOM     17  C   *   A 117      41.790  54.079  33.432  1.00  0.00      A    C  
ATOM     18  N   *   A 118      41.594  54.934  32.652  1.00  0.00      A    N  
ATOM     19  H   *   A 119      43.694  52.248  32.951  1.00  0.00      A    H  
ATOM     20  H   *   A 120      44.333  50.369  34.369  1.00  0.00      A    H  
ATOM     21  H   *   A 121      44.108  49.790  37.148  1.00  0.00      A    H  
ATOM     22  C   *   A 122      42.146  49.054  37.737  1.00  0.00      A    C  
ATOM     23  C   *   A 123      42.675  48.761  39.003  1.00  0.00      A    C  
ATOM     24  C   *   A 124      41.859  48.278  39.998  1.00  0.00      A    C """
        )


def test_process_pdb_file():
    # Define the input and output file paths
    original_cwd = Path(os.getcwd())
    input_pdb_filename = test_data_directory / "0_unk_hoh.pdb"

    shutil.copy(str(input_pdb_filename), ".")

    # Process the provided PDB file
    process_pdb_file(input_pdb_filename)

    # Read the modified output PDB file
    with open(input_pdb_filename, "r") as f:
        modified_data = f.read()

    # Check if the modified data contains the expected residues
    assert "HOH" in modified_data
    assert "UNK" in modified_data


def test_renumber_atoms_in_residues(sample_pdb_data, tmp_path):
    input_pdb_filename = tmp_path / "input.pdb"
    output_pdb_filename = tmp_path / "output.pdb"

    # Create a mock PDB file
    input_pdb_filename.write_text(
        """ATOM      1  N   UNK A 454      43.493  48.319  35.835  1.00  0.00      A    N  
ATOM      2  N1  UNK A 454      44.740  47.862  35.697  1.00  0.00      A    N  
ATOM      3  C14 UNK A 454      44.608  46.866  34.829  1.00  0.00      A    C  
ATOM      4  N2  UNK A 454      43.265  46.644  34.450  1.00  0.00      A    N  
ATOM      5  C7  UNK A 454      42.607  47.556  35.077  1.00  0.00      A    C  
ATOM      6  H5  UNK A 454      41.542  47.701  34.954  1.00  0.00      A    H  
ATOM      7  H10 UNK A 454      45.308  46.132  34.453  1.00  0.00      A    H  
ATOM      8  C   UNK A 454      43.168  49.513  36.656  1.00  0.00      A    C  
ATOM      9  C2  UNK A 454      42.743  50.705  35.818  1.00  0.00      A    C  
ATOM     10  C4  UNK A 454      43.545  51.052  34.671  1.00  0.00      A    C  
ATOM     11  C9  UNK A 454      43.171  52.151  33.897  1.00  0.00      A    C  
ATOM     12  C13 UNK A 454      42.090  52.924  34.222  1.00  0.00      A    C  
ATOM     13  C11 UNK A 454      41.393  52.671  35.378  1.00  0.00      A    C  
ATOM     14  C6  UNK A 454      41.793  51.635  36.268  1.00  0.00      A    C  
ATOM     15  H4  UNK A 454      41.220  51.358  37.148  1.00  0.00      A    H  
ATOM     16  H9  UNK A 454      40.518  53.291  35.552  1.00  0.00      A    H  
ATOM     17  C16 UNK A 454      41.790  54.079  33.432  1.00  0.00      A    C  
ATOM     18  N4  UNK A 454      41.594  54.934  32.652  1.00  0.00      A    N  
ATOM     19  H7  UNK A 454      43.694  52.248  32.951  1.00  0.00      A    H  
ATOM     20  H2  UNK A 454      44.333  50.369  34.369  1.00  0.00      A    H  
ATOM     21  H   UNK A 454      44.108  49.790  37.148  1.00  0.00      A    H  
ATOM     22  C1  UNK A 454      42.146  49.054  37.737  1.00  0.00      A    C  
ATOM     23  C5  UNK A 454      42.675  48.761  39.003  1.00  0.00      A    C  
ATOM     24  C10 UNK A 454      41.859  48.278  39.998  1.00  0.00      A    C
ATOM     25  H8  UNK A 454      42.284  48.099  40.981  1.00  0.00      A    H  
ATOM     26  H3  UNK A 454      43.752  48.806  39.135  1.00  0.00      A    H  
ATOM     27  C3  UNK A 454      40.774  48.885  37.463  1.00  0.00      A    C  
ATOM     28  H1  UNK A 454      40.310  49.079  36.500  1.00  0.00      A    H  
ATOM     29  C8  UNK A 454      39.907  48.435  38.509  1.00  0.00      A    C  
ATOM     30  H6  UNK A 454      38.833  48.310  38.406  1.00  0.00      A    H  
ATOM     31  C12 UNK A 454      40.466  48.125  39.823  1.00  0.00      A    C  
ATOM     32  C15 UNK A 454      39.627  47.605  40.833  1.00  0.00      A    C  
ATOM     33  N3  UNK A 454      38.981  47.235  41.740  1.00  0.00      A    N"""
    )

    renumber_atoms_in_residues(str(input_pdb_filename), str(output_pdb_filename), "UNK")
    assert output_pdb_filename.exists()


@pytest.fixture
def sample_pdb_info():
    return """
ATOM   741  N   UNK A 454      43.056  48.258  36.260  1.00  0.00      LIG  X  
ATOM   742  N1  UNK A 454      44.324  47.906  35.996  1.00  0.00      LIG  X  
ATOM   743  C14 UNK A 454      44.132  46.990  35.061  1.00  0.00      LIG  X  
    """


def test_process_pdb(sample_pdb_info):
    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as temp_file:
        temp_filename = temp_file.name
        temp_file.write(sample_pdb_info)

    print("Temp Data:")
    print(temp_filename)
    output_filename = "output_pdb_test.pdb"
    process_pdb(temp_filename, output_filename)

    with open(output_filename, "r") as f:
        modified_data = f.read()

    print("Modified Data:")
    print(modified_data)

    assert " LIG  N" in modified_data
    assert " LIG  C" in modified_data
    assert " LIG  X" not in modified_data

    # Clean up temporary and output files
    os.remove(temp_filename)
    os.remove(output_filename)


def test_extract_and_save_ligand_as_sdf():
    input_pdb_filename = topology_metal
    output_filename = "ligand_changed.sdf"
    target_resname = ligand_resname

    extract_and_save_ligand_as_sdf(input_pdb_filename, output_filename, target_resname)

    assert output_filename is not None
    os.remove("ligand_changed.sdf")


test_data_directory = Path("openmmdl/tests/data/in")
TEST_LIGAND_FILE = f"{test_data_directory}/CVV.sdf"
TEST_OUTPUT_FILE = "CVV.smi"


def test_increase_ring_indices():
    # Test case 1: Check if ring indices are correctly increased
    ring = [1, 2, 3]
    lig_index = 10
    result = increase_ring_indices(ring, lig_index)
    assert result == [11, 12, 13]

    # Test case 2: Check with a different lig_index
    ring = [3, 4, 5]
    lig_index = 20
    result = increase_ring_indices(ring, lig_index)
    assert result == [23, 24, 25]


def test_convert_ligand_to_smiles():
    # Convert the ligand structure to SMILES in the same directory as the input SDF file
    convert_ligand_to_smiles(TEST_LIGAND_FILE, TEST_OUTPUT_FILE)

    # Verify that the output SMILES file was created in the same directory as the input file
    assert os.path.exists(TEST_OUTPUT_FILE)

    # Optionally, you can also read and validate the content of the output SMILES file
    with open(TEST_OUTPUT_FILE, "r") as smi_file:
        smiles_lines = smi_file.readlines()
        assert len(smiles_lines) > 0  # Check that there are SMILES representations
