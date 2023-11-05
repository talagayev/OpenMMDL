import numpy as np
import pandas as pd
import re
import os
import matplotlib.pyplot as plt
import pytest

from openmmdl.openmmdl_analysis.binding_mode_processing import *


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

# Define a sample test case
@pytest.fixture
def sample_data():
    # Create sample DataFrames for testing
    df = pd.DataFrame({'FRAME': [1, 2, 3, 4],
                       'Value1': [10, 20, 30, 40],
                       'Value2': [100, 200, 300, 400]})

    new_df = pd.DataFrame({'FRAME': [2, 3, 4],
                           'Value1': [21, 31, 41],
                           'Value2': [210, 310, 410]})

    unique_data = {'Value1': 'Value1', 'Value2': 'Value2'}

    return df, new_df, unique_data

def test_update_values(sample_data):
    df, new_df, unique_data = sample_data

    # Call the function to be tested
    update_values(df, new_df, unique_data)

    # Check if the values have been updated as expected
    assert df.loc[df['FRAME'] == 2, 'Value1'].values[0] == 21
    assert df.loc[df['FRAME'] == 2, 'Value2'].values[0] == 210
    assert df.loc[df['FRAME'] == 3, 'Value1'].values[0] == 31
    assert df.loc[df['FRAME'] == 3, 'Value2'].values[0] == 310
    assert df.loc[df['FRAME'] == 4, 'Value1'].values[0] == 41
    assert df.loc[df['FRAME'] == 4, 'Value2'].values[0] == 410

# Define a test case that uses the fixture
def test_remove_duplicate_values(test_remove_duplicates_data):
    input_data, expected_output = test_remove_duplicates_data
    assert remove_duplicate_values(input_data) == expected_output
    
