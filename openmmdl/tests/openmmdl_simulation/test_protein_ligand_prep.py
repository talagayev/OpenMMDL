import pytest
import os
from pathlib import Path
from openmmdl.openmmdl_simulation.scripts.protein_ligand_prep import *

# Print current working directory
print("Current working directory:", os.getcwd())

# Define the full path to the input file
input_pdb_filename = "openmmdl/tests/data/in/0_unk_hoh.pdb"
print("Full path to input file:", os.path.abspath(input_pdb_filename))

test_data_directory = Path("openmmdl/tests/data/in")

# Define full paths to test files
TEST_LIGAND_FILE = test_data_directory / 'lig.sdf'
TEST_PROTEIN = test_data_directory / '0_unk_hoh.pdb'

# Test the protein_choice function
def test_protein_choice():
    prepared_protein = protein_choice("Yes", TEST_PROTEIN)
    assert isinstance(prepared_protein, pdbfixer.PDBFixer)

# Test the prepare_ligand function
def test_prepare_ligand():
    rdkitmolh = prepare_ligand(TEST_LIGAND_FILE)
    assert isinstance(rdkitmolh, Chem.Mol)

# Test the rdkit_to_openmm function
def test_rdkit_to_openmm():
    # Replace 'valid_smiles_string' with an actual valid SMILES string
    valid_smiles_string = '...'
    rdkit_mol = Chem.MolFromSmiles(valid_smiles_string)
    omm_molecule = rdkit_to_openmm(rdkit_mol, 'Ligand')
    assert isinstance(omm_molecule, app.Modeller)

# Other test functions for the remaining functions

if __name__ == '__main__':
    pytest.main()
