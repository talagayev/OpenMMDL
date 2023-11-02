import pytest
import os
import rdkit
from rdkit import Chem
from simtk import openmm, unit
from pathlib import Path
import openmm
from openmm.app import PDBFile
from pdbfixer import PDBFixer
from openmmdl.openmmdl_simulation.scripts.protein_ligand_prep import *

# Print current working directory
print("Current working directory:", os.getcwd())

# Assuming that 'test_data_directory' is properly defined in your test setup
test_data_directory = "openmmdl/tests/data/in"


test_data_directory = Path("openmmdl/tests/data/in")


# Define the full path to the input SDF file
TEST_LIGAND_FILE = f"{test_data_directory}/CVV.sdf"
TEST_MOL_FILE = f"{test_data_directory}/CVV.mol"
TEST_MOL2_FILE = f"{test_data_directory}/CVV.mol2"
TEST_PROTEIN = test_data_directory / '6b73.pdb'

# Test the protein_choice function
def test_protein_choice():
    prepared_protein = protein_choice("Yes", TEST_PROTEIN)
    assert isinstance(prepared_protein, pdbfixer.PDBFixer)

# Test the prepare_ligand function
def test_prepare_ligand():
    # Test the function with the sample ligand file.
    rdkit_mol_sdf = prepare_ligand(TEST_LIGAND_FILE, minimize_molecule=True)
    rdkit_mol_mol = prepare_ligand(TEST_MOL_FILE, minimize_molecule=False)
    rdkit_mol_mol2 = prepare_ligand(TEST_MOL2_FILE, minimize_molecule=False)
    
    # Add your assertions here to check if the preparation worked as expected
    assert rdkit_mol_sdf is not None  # Check if the result is not None
    assert rdkit_mol_mol is not None  # Check if the result is not None
    assert rdkit_mol_mol2 is not None  # Check if the result is not None



def test_rdkit_to_openmm_conversion():
    ligand_prepared = prepare_ligand(TEST_LIGAND_FILE, minimize_molecule=False)
    
    # Convert the RDKit molecule to an OpenMM modeller object
    omm_ligand = rdkit_to_openmm(ligand_prepared, 'UNK')

    # Check if the OpenMM modeller is an instance of OpenMM's app.Modeller
    assert omm_ligand is not None

def test_water_conversion():
    # Load the sample PDB file
    protein_name = "6b73.pdb"
    
    # Create a sample PDBFixer object and specify a model water
    fixer = PDBFixer()
    model_water = "TIP4P-EW"

    # Call the water_conversion function
    modeller = water_conversion(model_water, fixer, protein_name)

    # Check if the modeller object is an instance of Modeller
    assert isinstance(modeller, Modeller)

if __name__ == '__main__':
    pytest.main()
