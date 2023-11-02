import pytest
import os
from pathlib import Path
import mdtraj as md

from openmmdl.openmmdl_simulation.scripts.post_md_conversions import mdtraj_conversion, MDanalysis_conversion



test_data_directory = Path("openmmdl/tests/data/in")
pdb_file = f"{test_data_directory}/0_unk_hoh.pdb"

def test_mdtraj_conversion():
    # Create temporary directories to save the output files
    output_file_dcd = f"{test_data_directory}/0_unk_hoh.pdb"
    output_file_xtc = 'centered_old_coordinates.xtc'
    output_file_pdb = 'centered_old_coordinates_top.pdb'
    output_file_gro = 'centered_old_coordinates_top.gro'
  
    mdtraj_conversion(pdb_file, "gro_xtc")
    mdtraj_conversion(pdb_file, "pdb_dcd")
    
    assert output_file_dcd is not None
    assert output_file_xtc is not None
    assert output_file_pdb is not None
    assert output_file_gro is not None

