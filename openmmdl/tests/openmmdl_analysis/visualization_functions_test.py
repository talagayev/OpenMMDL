import numpy as np
import pandas as pd
import re
import shutil
import subprocess
import os
import matplotlib.pyplot as plt
import pytest
from openmmdl.openmmdl_analysis.visualization_functions import *



# visualization_functions tests
@pytest.fixture
def sample_dataframe_interacting_water_ids():
    data = {
        'Interaction1': [0, 1, 0, 1, 0],
        'Interaction2': [1, 0, 0, 0, 1],
        'WATER_IDX': [101, 102, None, 104, 105],
        'FRAME': [1, 2, 3, 4, 5]  
    }
    df_all = pd.DataFrame(data)
    return df_all

def test_interacting_water_ids(sample_dataframe_interacting_water_ids):
    waterbridge_interactions = ['Interaction1', 'Interaction2']
    
    result = interacting_water_ids(sample_dataframe_interacting_water_ids, waterbridge_interactions)

    expected_interacting_waters = [101, 102, 104, 105]

    assert sorted(result) == sorted(expected_interacting_waters)
    

@pytest.fixture
def sample_dataframe_cloud_json_generation():
    data = {
        'LIGCOO': [
            "(1.0, 2.0, 3.0)",
            "(4.0, 5.0, 6.0)",
            "(7.0, 8.0, 9.0)",
        ],
        'INTERACTION': [
            'hydrophobic',
            'acceptor',
            'donor',
        ],
        'PROTISDON': [
            'False',
            'True',
            'False',
        ],
        'PROTISPOS': [
            'False',
            'False',
            'True',
        ],
    }
    df_all = pd.DataFrame(data)
    return df_all

#def test_cloud_json_generation(sample_dataframe_cloud_json_generation):
#    result = cloud_json_generation(sample_dataframe_cloud_json_generation)
#
#    expected_clouds = {
#        'acceptor': {
#            "coordinates": [[4.0, 5.0, 6.0]],
#            "color": [1.0, 0.0, 0.0],
#            "radius": 0.05,
#        },
#        'donor': {
#            "coordinates": [[7.0, 8.0, 9.0]],
#            "color": [0.0, 1.0, 0.0],
#            "radius": 0.05,
#        },
#        'hydrophobic': {
#            "coordinates": [[1.0, 2.0, 3.0]],
#            "color": [1.0, 1.0, 0.0],
#            "radius": 0.05,
#        },
#        'waterbridge': {
#            "coordinates": [],
#            "color": [0.0, 1.0, 0.9],
#            "radius": 0.05,
#        },
#        'negative_ionizable': {
#            "coordinates": [],
#            "color": [0.0, 0.0, 1.0],
#            "radius": 0.05,
#        },
#        'positive_ionizable': {
#            "coordinates": [],
#            "color": [1.0, 0.0, 0.0],
#            "radius": 0.05,
#        },
#        'pistacking': {
#            "coordinates": [],
#            "color": [0.0, 0.0, 1.0],
#            "radius": 0.05,
#        },
#        'pication': {
#            "coordinates": [],
#            "color": [0.0, 0.0, 1.0],
#            "radius": 0.05,
#        },
#        'halogen': {
#            "coordinates": [],
#            "color": [1.0, 0.0, 0.9],
#            "radius": 0.05,
#        },
#        'metal': {
#            "coordinates": [],
#            "color": [1.0, 0.6, 0.0],
#            "radius": 0.05,
#        },
#    }
#
#    assert result == expected_clouds

def test_run_visualization(tmp_path, monkeypatch):
    # tmp_path is a pytest fixture that provides a temporary directory
    
    # Set up the temporary directory
    tmp_dir = tmp_path / "test_dir"
    tmp_dir.mkdir()
    
    # Create a mock visualization notebook for testing
    mock_notebook_path = tmp_dir / 'visualization.ipynb'
    mock_notebook_content = "Mock notebook content"
    with open(mock_notebook_path, 'w') as mock_notebook:
        mock_notebook.write(mock_notebook_content)
    
    # Mock the os.path.dirname function to return the temporary directory
    monkeypatch.setattr(os.path, 'dirname', lambda path: str(tmp_dir))
    
    # Mock the os.getcwd function to return the temporary directory
    monkeypatch.setattr(os, 'getcwd', lambda: str(tmp_dir))
    
    # Mock the subprocess.run function to check if it is called with the correct arguments
    with pytest.raises(subprocess.CalledProcessError):
        run_visualization()

    # Check if the notebook was copied to the current directory
    copied_notebook_path = tmp_dir / 'visualization.ipynb'
    assert copied_notebook_path.is_file()
    
    # Check if the content of the copied notebook is the same as the mock notebook
    with open(copied_notebook_path, 'r') as copied_notebook:
        assert copied_notebook.read() == mock_notebook_content
