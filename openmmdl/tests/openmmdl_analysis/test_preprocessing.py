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
def temp_pdb_file(tmp_path):
    input_pdb_filename = tmp_path / "test_input.pdb"
    # Copy the content of the provided PDB file to the temporary test file
    with open(pdb_file_path, "r") as src_pdb, open(input_pdb_filename, "w") as dest_pdb:
        dest_pdb.write(src_pdb.read())
    return input_pdb_filename

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

def test_convert_pdb_to_sdf(tmp_path):
    input_pdb_filename = tmp_path / "input.pdb"
    output_sdf_filename = tmp_path / "output.sdf"
    
    # Create a mock PDB file
    input_pdb_filename.write_text("""ATOM      1  N   UNK A 454      43.493  48.319  35.835  1.00  0.00      A    N  
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
ATOM     33  N3  UNK A 454      38.981  47.235  41.740  1.00  0.00      A    N""")

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
    # Create a mock PDB file with the updated input data
    input_data = """ATOM      1  N   UNK A 454      43.493  48.319  35.835  1.00  0.00      A    N  
ATOM      2  N1  LIG A 454      44.740  47.862  35.697  1.00  0.00      LIG  X  
ATOM      3  C14 LIG A 454      44.608  46.866  34.829  1.00  0.00      LIG  X  
ATOM      4  N2  LIG A 454      43.265  46.644  34.450  1.00  0.00      LIG  X  
ATOM      5  C7  LIG A 454      42.607  47.556  35.077  1.00  0.00      LIG  X  
ATOM      6  H5  LIG A 454      41.542  47.701  34.954  1.00  0.00      LIG  X  
ATOM      7  H10 LIG A 454      45.308  46.132  34.453  1.00  0.00      LIG  X  
ATOM      8  C   LIG A 454      43.168  49.513  36.656  1.00  0.00      LIG  X  
ATOM      9  C2  LIG A 454      42.743  50.705  35.818  1.00  0.00      LIG  X  
ATOM     10  C4  LIG A 454      43.545  51.052  34.671  1.00  0.00      LIG  X  
ATOM     11  C9  LIG A 454      43.171  52.151  33.897  1.00  0.00      LIG  X  
ATOM     12  C13 LIG A 454      42.090  52.924  34.222  1.00  0.00      LIG  X  
ATOM     13  C11 LIG A 454      41.393  52.671  35.378  1.00  0.00      LIG  X  
ATOM     14  C6  LIG A 454      41.793  51.635  36.268  1.00  0.00      LIG  X  
ATOM     15  H4  LIG A 454      41.220  51.358  37.148  1.00  0.00      LIG  X  
ATOM     16  H9  LIG A 454      40.518  53.291  35.552  1.00  0.00      LIG  X  
ATOM     17  C16 LIG A 454      41.790  54.079  33.432  1.00  0.00      LIG  X  
ATOM     18  N4  LIG A 454      41.594  54.934  32.652  1.00  0.00      LIG  X  
ATOM     19  H7  LIG A 454      43.694  52.248  32.951  1.00  0.00      LIG  X  
ATOM     20  H2  LIG A 454      44.333  50.369  34.369  1.00  0.00      LIG  X  
ATOM     21  H   LIG A 454      44.108  49.790  37.148  1.00  0.00      LIG  X  
ATOM     22  C1  LIG A 454      42.146  49.054  37.737  1.00  0.00      LIG  X  
ATOM     23  C5  LIG A 454      42.675  48.761  39.003  1.00  0.00      LIG  X  
ATOM     24  C10 LIG A 454      41.859  48.278  39.998  1.00  0.00      LIG  X
ATOM     25  H8  LIG A 454      42.284  48.099  40.981  1.00  0.00      LIG  X  
ATOM     26  H3  LIG A 454      43.752  48.806  39.135  1.00  0.00      LIG  X  
ATOM     27  C3  LIG A 454      40.774  48.885  37.463  1.00  0.00      LIG  X  
ATOM     28  H1  LIG A 454      40.310  49.079  36.500  1.00  0.00      LIG  X  
ATOM     29  C8  LIG A 454      39.907  48.435  38.509  1.00  0.00      LIG  X  
ATOM     30  H6  LIG A 454      38.833  48.310  38.406  1.00  0.00      LIG  X  
ATOM     31  C12 LIG A 454      40.466  48.125  39.823  1.00  0.00      LIG  X  
ATOM     32  C15 LIG A 454      39.627  47.605  40.833  1.00  0.00      LIG  X  
ATOM     33  N3  LIG A 454      38.981  47.235  41.740  1.00  0.00      LIG  X """

    # Write the input data to the mock input file
    input_file.write_text(input_data)

    # Apply the replace_atom_type function to the input data
    modified_data = replace_atom_type(input_data)

    # Write the modified data to the mock output file
    output_file.write_text(modified_data)

    # Check if the modification was successful
    assert 'LIG  C' in output_file.read_text()
