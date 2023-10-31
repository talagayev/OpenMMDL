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
    
    assert Chem.MolToSmiles(rdkit_mol)  # Check if a valid SMILES can be generated.

    # Check if hydrogen atoms are added.
    assert rdkit_mol.GetNumAtoms() > 0
    assert rdkit_mol.GetNumAtoms() > rdkit_mol.GetNumHeavyAtoms()

    # Check if chiral tags are assigned.
    assert all(atom.HasChiralTag() for atom in rdkit_mol.GetAtoms())

    # Check if minimization was performed when selected.
    assert isinstance(rdkit_mol, Chem.Mol)  # Check if the molecule is still an RDKit molecule.

    # Check the conversion to an OpenFF Molecule object.
    openff_mol = Molecule(rdkit_mol)
    assert isinstance(openff_mol, Molecule)

def test_rdkit_to_openmm():
    # Create an RDKit molecule (rdkit_mol) and provide a name.
    rdkit_mol = Chem.MolFromSmiles('CCO')  # Example molecule: Ethanol
    name = "Ethanol"

    # Convert the RDKit molecule to an OpenMM Modeller object.
    omm_mol = rdkit_to_openmm(rdkit_mol, name)

    # Check if the returned object is an instance of openmm.app.Modeller.
    assert isinstance(omm_mol, openmm.app.Modeller)

    # Check if the Modeller object's topology contains the expected number of atoms.
    num_atoms_expected = rdkit_mol.GetNumAtoms()
    num_atoms_actual = omm_mol.topology.getNumAtoms()
    assert num_atoms_expected == num_atoms_actual

    # Check if the Modeller object's name matches the provided name.
    assert omm_mol.topology.getPeriodicBoxVectors() == name

    # Check if the Modeller object's positions are in nanometers.
    positions = omm_mol.getPositions()
    for position in positions:
        assert position.unit == unit.nanometers

if __name__ == '__main__':
    pytest.main()
