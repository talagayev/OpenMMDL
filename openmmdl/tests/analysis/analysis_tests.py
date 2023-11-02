import numpy as np
import pandas as pd
import re
import os
import matplotlib.pyplot as plt
import pytest
from openmmdl.openmmdl_analysis.barcode_generation import *
from openmmdl.openmmdl_analysis.binding_mode_processing import *
from openmmdl.openmmdl_analysis.pml_writer import *
from openmmdl.openmmdl_analysis.visualization_functions import *


# Barcode generation tests
@pytest.fixture
def sample_dataframe_barcode_generation():
    data = {
        'FRAME': [1, 1, 2, 2, 3],
        'Interaction1': [1, 0, 1, 0, 0],
        'Interaction2': [0, 0, 0, 1, 1],
        'WATER_IDX': [101, 102, 103, 104, 105],
    }
    return pd.DataFrame(data)

def test_barcodegeneration(sample_dataframe_barcode_generation):
    interaction = 'Interaction1'
    barcode = barcodegeneration(sample_dataframe_barcode_generation, interaction)
    
    assert isinstance(barcode, np.ndarray)
    
    expected_barcode = np.array([1, 1, 0])
    assert np.array_equal(barcode, expected_barcode)
    
def test_waterids_barcode_generator(sample_dataframe_barcode_generation):
    interaction = 'Interaction2'
    waterid_barcode = waterids_barcode_generator(sample_dataframe_barcode_generation, interaction)
    
    # Test if the output is a list
    assert isinstance(waterid_barcode, list)
    
    # Test the expected waterid barcode for the sample dataframe and interaction
    expected_waterid_barcode = [0, 104, 105]
    assert waterid_barcode == expected_waterid_barcode


# binding_mode_processing tests
@pytest.fixture
def sample_dataframe_bindingmode_processing():
    data = {
        'FRAME': {0: 1, 1: 2, 2: 3, 3: 2},
        'Prot_partner': {0: 'A', 1: 'B', 2: 'C', 3: 'A'},
        'INTERACTION': {0: 'hydrophobic', 1: 'hbond', 2: 'saltbridge', 3: 'hydrophobic'},
        'LIGCARBONIDX': {0: 101, 1: 102, 2: 103, 3: 102},
        'DONORIDX': {0: 201, 1: 202, 2: 203, 3: 202},
        'ACCEPTORIDX': {0: 301, 1: 302, 2: 303, 3: 302},
        'PROTISDON': {0: True, 1: False, 2: True, 3: False},
        'LIG_IDX_LIST': {0: [1, 2], 1: [3, 4], 2: [5, 6], 3: [3, 4]},
        'LIG_GROUP': {0: 'Group1', 1: 'Group2', 2: 'Group3', 3: 'Group1'},
        'PROTISPOS': {0: True, 1: False, 2: True, 3: True}
    }
    return pd.DataFrame(data)


def test_gather_interactions(sample_dataframe_bindingmode_processing):
    df = sample_dataframe_bindingmode_processing
    ligand_rings = [[101], [102], [103]]  # Define sample ligand rings for testing

    result = gather_interactions(df, ligand_rings)

    # Assert that the result is a dictionary
    assert isinstance(result, dict)

    # Check specific values in the generated dictionary for known interactions based on the updated fixture
    expected_result = {
        1: {0: 'A_101_hydrophobic'},
        2: {1: 'B_202_Donor_hbond', 3: 'A_102_hydrophobic'},
        3: {2: 'C_[5, 6]_Group3_NI_saltbridge'}
    }

    # Check if the actual result matches the expected result
    assert result == expected_result

@pytest.fixture
def test_remove_duplicates_data():
    input_data = {
        'a': {'x': 1, 'y': 2, 'z': 1},
        'b': {'p': 3, 'q': 3, 'r': 4}
    }
    expected_output = {
        'a': {'x': 1, 'y': 2},
        'b': {'p': 3, 'r': 4}
    }
    return input_data, expected_output

# Define a test case that uses the fixture
def test_remove_duplicate_values(test_remove_duplicates_data):
    input_data, expected_output = test_remove_duplicates_data
    assert remove_duplicate_values(input_data) == expected_output
    

# pml_writer tests
@pytest.fixture
def sample_dataframe_generate_pharmacophore_centers():
    data = {
        'Hydrophobic': [1, 1, 0, 1, 0],
        'Ionic': [0, 1, 0, 0, 1],
        'LIGCOO': ["(1.0, 2.0, 3.0)", "(2.0, 3.0, 4.0)", "(3.0, 4.0, 5.0)", "(4.0, 5.0, 6.0)", "(5.0, 6.0, 7.0)"]
    }
    df = pd.DataFrame(data)
    return df

@pytest.fixture
def sample_interactions_generate_pharmacophore_centers():
    return ['Hydrophobic', 'Ionic']

def test_generate_pharmacophore_centers(sample_dataframe_generate_pharmacophore_centers, sample_interactions_generate_pharmacophore_centers):
    result = generate_pharmacophore_centers(sample_dataframe_generate_pharmacophore_centers, sample_interactions_generate_pharmacophore_centers)
    
    expected_pharmacophore = {
        'Hydrophobic': [2.333, 3.333, 4.333],
        'Ionic': [3.5, 4.5, 5.5]
    }
    
    assert result == expected_pharmacophore


@pytest.fixture
def sample_dataframe_generate_pharmacophore_vectors():
    # Create a sample dataframe for testing
    data = {
        'HBDonors': [1, 0, 1, 0, 1],
        'HBAcceptors': [0, 1, 0, 1, 0],
        'LIGCOO': [
            "(1.0, 2.0, 3.0)",
            "(2.0, 3.0, 4.0)",
            "(3.0, 4.0, 5.0)",
            "(4.0, 5.0, 6.0)",
            "(5.0, 6.0, 7.0)"
        ],
        'PROTCOO': [
            "(0.5, 1.5, 2.5)",
            "(1.5, 2.5, 3.5)",
            "(2.5, 3.5, 4.5)",
            "(3.5, 4.5, 5.5)",
            "(4.5, 5.5, 6.5)"
        ]
    }
    df = pd.DataFrame(data)
    return df

@pytest.fixture
def sample_interactions_generate_pharmacophore_vectors():
    return ['HBDonors', 'HBAcceptors']

def test_generate_pharmacophore_vectors(sample_dataframe_generate_pharmacophore_vectors, sample_interactions_generate_pharmacophore_vectors):
    result = generate_pharmacophore_vectors(sample_dataframe_generate_pharmacophore_vectors, sample_interactions_generate_pharmacophore_vectors)
    
    expected_pharmacophore = {
        'HBDonors': [
            [3.0, 4.0, 5.0],
            [2.5, 3.5, 4.5]
        ],
        'HBAcceptors': [
            [3.0, 4.0, 5.0],
            [2.5, 3.5, 4.5]
        ]
    }

    assert result == expected_pharmacophore
    

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

def test_cloud_json_generation(sample_dataframe_cloud_json_generation):
    result = cloud_json_generation(sample_dataframe_cloud_json_generation)

    expected_clouds = {
        'acceptor': {
            "coordinates": [[4.0, 5.0, 6.0]],
            "color": [1.0, 0.0, 0.0],
            "radius": 0.05,
        },
        'donor': {
            "coordinates": [[7.0, 8.0, 9.0]],
            "color": [0.0, 1.0, 0.0],
            "radius": 0.05,
        },
        'hydrophobic': {
            "coordinates": [[1.0, 2.0, 3.0]],
            "color": [1.0, 1.0, 0.0],
            "radius": 0.05,
        },
        'waterbridge': {
            "coordinates": [],
            "color": [0.0, 1.0, 0.9],
            "radius": 0.05,
        },
        'negative_ionizable': {
            "coordinates": [],
            "color": [0.0, 0.0, 1.0],
            "radius": 0.05,
        },
        'positive_ionizable': {
            "coordinates": [],
            "color": [1.0, 0.0, 0.0],
            "radius": 0.05,
        },
        'pistacking': {
            "coordinates": [],
            "color": [0.0, 0.0, 1.0],
            "radius": 0.05,
        },
        'pication': {
            "coordinates": [],
            "color": [0.0, 0.0, 1.0],
            "radius": 0.05,
        },
        'halogen': {
            "coordinates": [],
            "color": [1.0, 0.0, 0.9],
            "radius": 0.05,
        },
        'metal': {
            "coordinates": [],
            "color": [1.0, 0.6, 0.0],
            "radius": 0.05,
        },
    }

    assert result == expected_clouds