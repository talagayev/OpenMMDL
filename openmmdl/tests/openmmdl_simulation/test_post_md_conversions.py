import pytest
import os
import mdtraj as md

from openmmdl.openmmdl_simulation.scripts.post_md_conversions import mdtraj_conversion, MDanalysis_conversion


pdb_file = "0_unk_hoh.pdb"
test_data_dir = "openmmdl/tests/data/in"


def test_mdtraj_conversion():
    # Create temporary directories to save the output files
    output_file_dcd = 'centered_old_coordinates.dcd'
    output_file_xtc = 'centered_old_coordinates.xtc'
    output_file_pdb = 'centered_old_coordinates_top.pdb'
    output_file_gro = 'centered_old_coordinates_top.gro'
    
    original_cwd = os.getcwd()
    os.chdir(test_data_dir)
  
    mdtraj_conversion(pdb_file, "gro_xtc")
    mdtraj_conversion(pdb_file, "pdb_dcd")
    
    assert output_file_dcd is not None
    assert output_file_xtc is not None
    assert output_file_pdb is not None
    assert output_file_gro is not None
