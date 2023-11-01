import pytest
import os
import rdkit
from rdkit import Chem
from simtk import openmm, unit
from pathlib import Path
from openmmdl.openmmdl_simulation.scripts.protein_ligand_prep import *

# Print current working directory
print("Current working directory:", os.getcwd())

# Assuming that 'test_data_directory' is properly defined in your test setup
test_data_directory = "openmmdl/tests/data/in"


test_data_directory = Path("openmmdl/tests/data/in")


# Define the full path to the input SDF file
TEST_LIGAND_FILE = f"{test_data_directory}/CVV.sdf"
TEST_PROTEIN = test_data_directory / '6b73.pdb'

# Test the protein_choice function
def test_protein_choice():
    prepared_protein = protein_choice("Yes", TEST_PROTEIN)
    assert isinstance(prepared_protein, pdbfixer.PDBFixer)

# Test the prepare_ligand function
def test_prepare_ligand():
    # Test the function with the sample ligand file.
    rdkit_mol = prepare_ligand(TEST_LIGAND_FILE, minimize_molecule=True)
    
    # Add your assertions here to check if the preparation worked as expected
    assert rdkit_mol is not None  # Check if the result is not None

# Define test cases
@pytest.mark.parametrize("rdkit_mol, name, expected_num_atoms", [
    (Chem.MolFromSmiles("CCO"), "Ethanol", 6),  # Example with ethanol molecule
    # Add more test cases as needed
])

def test_rdkit_to_openmm(rdkit_mol, name, expected_num_atoms):
    # Call the conversion function
    omm_molecule = rdkit_to_openmm(rdkit_mol, name)
    
    # Check if the result is an OpenMM Modeller
    assert isinstance(omm_molecule, app.Modeller)

    # Check the number of atoms in the OpenMM molecule
    assert len(omm_molecule.topology.atoms()) == expected_num_atoms

if __name__ == '__main__':
    pytest.main()
