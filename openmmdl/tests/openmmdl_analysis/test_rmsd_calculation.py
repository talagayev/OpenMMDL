import os
import pytest
from pathlib import Path
import pandas as pd
import numpy as np
import mdtraj as md

from openmmdl.openmmdl_analysis.rmsd_calculation import rmsd_for_atomgroups

test_data_directory = Path("openmmdl/tests/data/in")
topology_file = "0_unk_hoh.pdb"
trajectory_file = "all_50.dcd"
selection1 = "protein"
selection2 = ["ligand"]

def test_rmsd_for_atomgroups():

    # Go into folder with data
    original_cwd = os.getcwd()
    os.chdir(test_data_directory)

    # Call the function
    rmsd_df = rmsd_for_atomgroups(topology_file, trajectory_file, selection1, selection2)

    # Check if the output DataFrame has the correct structure
    assert isinstance(rmsd_df, pd.DataFrame)
    assert rmsd_df.index.name == "frame"
    
    # Check if the CSV file exists
    assert os.path.exists("RMSD_over_time.csv")
    
    # Check if the plot file exists
    assert os.path.exists("RMSD_over_time.png")
    
    # Cleanup created files after the test
    os.remove("RMSD_over_time.csv")
    os.remove("RMSD_over_time.png")
    os.chdir(original_cwd)
