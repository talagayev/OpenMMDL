import os
import pytest
from pathlib import Path
import pandas as pd
import numpy as np
import mdtraj as md
from plip.structure.preparation import PDBComplex, LigandFinder, Mol, PLInteraction

from openmmdl.openmmdl_analysis.interaction_gathering import characterize_complex, retrieve_plip_interactions


test_data_directory = Path("openmmdl/tests/data/in")
topology_file = f"{test_data_directory}/complex.pdb"
binding_site_id = "UNK:X:0"
lig_name = "UNK"

# Test the function
def test_characterize_complex():

    # Call the function
    interaction_set = characterize_complex(topology_file, binding_site_id)

    # Check if the function returns a PLInteraction object
    assert isinstance(interaction_set, PLInteraction)

def test_retrieve_plip_interactions(prepare_test_data):
    pdb_file, ligand_name = prepare_test_data

    # Call the function
    interactions = retrieve_plip_interactions(topology_file, lig_name)

    # Check if the function returns a dictionary
    assert isinstance(interactions, dict)
