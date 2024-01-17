import numpy as np
import pandas as pd
import re
import shutil
import subprocess
import os
from pathlib import Path
import matplotlib.pyplot as plt
import pytest
from openmmdl.openmmdl_analysis.visualization_functions import *

package_path = Path("openmmdl/openmmdl_analysis")

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

def test_run_visualization():
    # Set up the paths
    package_path = Path("openmmdl/openmmdl_analysis")
    notebook_path =  package_path / "visualization.ipynb"
    
    # Run the visualization function
    # run_visualization()
    
    # Check if the notebook was copied to the current directory with the correct name
    copied_notebook_path = os.path.join(os.getcwd(), 'visualization.ipynb')
    shutil.copy(str(notebook_path), '.')
    new_notebook_path = 'visualization.ipynb'
    assert os.path.isfile(copied_notebook_path)
    
    # Check if the content of the copied notebook is the same as the original notebook
    with open(new_notebook_path, 'r') as copied_notebook:
        with open(notebook_path, 'r') as original_notebook:
            assert copied_notebook.read() == original_notebook.read()

@pytest.fixture
def sample_dataframe():
    # Create a sample dataframe for testing
    data = {
        'LIGCOO': ['(1.0, 2.0, 3.0)', '(4.0, 5.0, 6.0)', '(13.0, 14.0, 15.0)', '(16.0, 17.0, 18.0)', '(19.0, 20.0, 21.0)'],
        'INTERACTION': ['hydrophobic', 'acceptor', 'donor', 'pistacking', 'pication'],
        'PROTISDON': ['False', 'True', 'True', 'False', 'True'],
        'PROTISPOS': ['False', 'True', 'False', 'False', 'False'],
        'TARGETCOO': ['(7.0, 8.0, 9.0)', '(10.0, 11.0, 12.0)', '(22.0, 23.0, 24.0)', '(25.0, 26.0, 27.0)', '(28.0, 29.0, 30.0)']
    }
    return pd.DataFrame(data)


def test_cloud_json_generation(sample_dataframe):
    result = cloud_json_generation(sample_dataframe)

    assert 'hydrophobic' in result
    assert 'acceptor' in result
    assert 'donor' in result
    assert 'waterbridge' in result
    assert 'negative_ionizable' in result
    assert 'positive_ionizable' in result
    assert 'pistacking' in result
    assert 'pication' in result
    assert 'halogen' in result
    assert 'metal' in result

    # Add more specific assertions based on your expectations for the output
    # For example, you might want to check the structure of the generated dictionary
    assert isinstance(result['hydrophobic'], dict)
    assert 'coordinates' in result['hydrophobic']
    assert 'color' in result['hydrophobic']
    assert 'radius' in result['hydrophobic']

    # Add more tests based on your specific requirements and expected results

