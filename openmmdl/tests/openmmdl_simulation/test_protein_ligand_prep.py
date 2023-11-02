import pytest
import os
import rdkit
from rdkit import Chem
import simtk.openmm.app as app
from simtk.openmm.app import PDBFile, Modeller
from simtk.openmm import unit
from simtk.openmm import Vec3
import mdtraj as md
import numpy as np
from pathlib import Path
import pdbfixer
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
    rdkit_mol_sdf = prepare_ligand(TEST_LIGAND_FILE, minimize_molecule=False)
    rdkit_mol_mol = prepare_ligand(TEST_MOL_FILE, minimize_molecule=False)
    rdkit_mol_mol2 = prepare_ligand(TEST_MOL2_FILE, minimize_molecule=False)
    
    # Add your assertions here to check if the preparation worked as expected
    assert rdkit_mol_sdf is not None  # Check if the result is not None
    assert rdkit_mol_mol is not None  # Check if the result is not None
    assert rdkit_mol_mol2 is not None  # Check if the result is not None


if __name__ == '__main__':
    pytest.main()
