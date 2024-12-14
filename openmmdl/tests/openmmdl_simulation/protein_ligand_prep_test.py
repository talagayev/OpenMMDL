import mdtraj as md
import numpy as np
import pdbfixer
from pdbfixer import PDBFixer
import simtk.openmm.app as app
import rdkit
from rdkit import Chem
from rdkit.Chem import rdForceFieldHelpers
from openff.toolkit.topology import Molecule
from simtk.openmm.app import PDBFile
from simtk.openmm import unit
from simtk.openmm import Vec3
import pytest

from openmmdl.openmmdl_simulation.protein_ligand_prep import PDBFormatter, MoleculePreparation
from openmmdl.openmmdl_simulation.parser import ConfigParser


# Fixture to create a test configuration file
@pytest.fixture
def create_config_file(tmp_path):
    def _create_config(content):
        config_file = tmp_path / "test_config.conf"
        config_file.write_text(content)
        return str(config_file)
    return _create_config

def test_pdbformatter(create_config_file):
    # Test case: No ligand
    config_file_no_ligand = create_config_file("""
protein = ../data/in/6b73.pdb
""")

    # Initialize the parser
    config_parser_no_ligand = ConfigParser(config_file_no_ligand)

    # Initialize the PDBFormatter
    formatter_no_ligand = PDBFormatter(config_parser_no_ligand)

    # Check that ligand is None and protein is located
    assert formatter_no_ligand.ligand is None
    assert formatter_no_ligand.protein == '../data/in/6b73.pdb'

    # Check the PDBFile case without ligand
    protein_pdb_no_ligand = formatter_no_ligand.get_protein_pdb()
    assert isinstance(protein_pdb_no_ligand, PDBFile)

    # Test case: With ligand
    config_file_with_ligand = create_config_file("""
ligand = ../data/in/lig.sdf
protein = ../data/in/6b73.pdb
""")

    # Initialize with ligand
    config_parser_with_ligand = ConfigParser(config_file_with_ligand)
    formatter_with_ligand = PDBFormatter(config_parser_with_ligand)

    # Check that the ligand was recognized
    assert formatter_with_ligand.ligand == '../data/in/lig.sdf'

    # 
    protein_pdb_with_ligand = formatter_with_ligand.get_protein_pdb()
    assert isinstance(protein_pdb_with_ligand, PDBFixer)

def test_prepare_ligand_rdkit_to_openmm():
    # Specify where the ligand is located
    ligand = '../data/in/lig.sdf'

    # Initialize ligand
    molecule_init = MoleculePreparation(ligand)

    # Prepare the ligand with RDKit
    molecule_prepared = molecule_init.prepare_ligand()
    assert isinstance(molecule_prepared, rdkit.Chem.rdchem.Mol)

    # Check that the ligand is converted to an Modeller object
    omm_lig = molecule_init.rdkit_to_openmm(molecule_prepared, 'UNK')
    assert isinstance(omm_lig, app.Modeller)
