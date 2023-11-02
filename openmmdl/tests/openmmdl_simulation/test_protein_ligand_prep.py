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
import openmm
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
    rdkit_mol_sdf = prepare_ligand(TEST_LIGAND_FILE, minimize_molecule=True)
    rdkit_mol_mol = prepare_ligand(TEST_MOL_FILE, minimize_molecule=False)
    rdkit_mol_mol2 = prepare_ligand(TEST_MOL2_FILE, minimize_molecule=False)
    
    # Add your assertions here to check if the preparation worked as expected
    assert rdkit_mol_sdf is not None  # Check if the result is not None
    assert rdkit_mol_mol is not None  # Check if the result is not None
    assert rdkit_mol_mol2 is not None  # Check if the result is not None



def test_rdkit_to_openmm_conversion():
    rdkit_mol = prepare_ligand(TEST_LIGAND_FILE, minimize_molecule=False)
    
    off_mol = Molecule.from_rdkit(rdkit_mol)

    # add name for molecule
    off_mol.name = name

    # add names for atoms
    element_counter_dict = {}
    for off_atom, rdkit_atom in zip(off_mol.atoms, rdkit_mol.GetAtoms()):
        element = rdkit_atom.GetSymbol()
        if element in element_counter_dict.keys():
            element_counter_dict[element] += 1
        else:
            element_counter_dict[element] = 1
        off_atom.name = element + str(element_counter_dict[element])
        
    # convert from OpenFF to OpenMM
    off_mol_topology = off_mol.to_topology()
    mol_topology = off_mol_topology.to_openmm()
    mol_positions = off_mol.conformers[0];
    new_mol_positions = []

    # convert units from Ångström to Nanometers
    for mol_position in off_mol.conformers[0]:
        new_mol_positions.append(mol_position.magnitude/10.0);
        print("yep")

    # combine topology and positions in modeller object
    omm_mol = app.Modeller(mol_topology, new_mol_positions * unit.nanometers)
    
    # Check if the OpenMM modeller is an instance of OpenMM's app.Modeller
    assert omm_mol is not None

def test_water_conversion():
    # Load the sample PDB file
    protein_name = "6b73.pdb"
    
    # Create a sample PDBFixer object and specify a model water
    modeller_pre_conversion = Modeller(PDBFile(str(TEST_PROTEIN)))  
    model_water = "TIP4P-EW"

    # Call the water_conversion function
    modeller = water_conversion(model_water, modeller_pre_conversion, protein_name)

    # Check if the modeller object is an instance of Modeller
    assert isinstance(modeller, Modeller)

if __name__ == '__main__':
    pytest.main()
