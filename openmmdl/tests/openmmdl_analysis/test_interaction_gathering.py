import os
import pytest
from pathlib import Path
import pandas as pd
import numpy as np
import mdtraj as md
from plip.structure.preparation import PDBComplex, LigandFinder, Mol, PLInteraction

from openmmdl.openmmdl_analysis.interaction_gathering import characterize_complex


test_data_directory = Path("openmmdl/tests/data/in")
topology_file = f"{test_data_directory}/0_unk_hoh.pdb"
binding_site_id = "UNK:X:0"

# Test the function
def test_characterize_complex():

    # Call the function
    interaction_set = characterize_complex(topology_file, binding_site_id)

    # Check if the function returns a PLInteraction object
    assert isinstance(interaction_set, PLInteraction)

    # Additional test cases and assertions can be added to test different scenarios.
    
    # Clean up created files after the test
    os.remove(pdb_file)
