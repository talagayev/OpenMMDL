import numpy as np
import pandas as pd
import rdkit
import MDAnalysis as mda
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
        'PROTISPOS': {0: True, 1: False, 2: True, 3: True},
        'DON_IDX': {0: 0, 1: 0, 2: 0, 3: 0},
        'DONORTYPE': {0: 0, 1: 0, 2: 0, 3: 0},
        'ACCEPTOR_IDX': {0: 0, 1: 0, 2: 0, 3: 0},
        'DONOR_IDX': {0: 0, 1: 0, 2: 0, 3: 0},
        'LOCATION': {0: 0, 1: 0, 2: 0, 3: 0},
        'METAL_IDX': {0: 0, 1: 0, 2: 0, 3: 0},
        'METAL_TYPE': {0: 0, 1: 0, 2: 0, 3: 0}
    }

    # Add 'halogen' and 'hbond' data to the existing DataFrame
    data['FRAME'][4] = 4  # Add a new 'FRAME' value
    data['Prot_partner'][4] = 'A'  # Add a new 'Prot_partner' value
    data['INTERACTION'][4] = 'halogen'  # Add 'halogen' interaction
    data['DON_IDX'][4] = 501  # DON_IDX for 'halogen'
    data['DONORTYPE'][4] = 'F'  # Halogen type
    data['ACCEPTOR_IDX'][4] = 0
    data['DONOR_IDX'][4] = 0
    data['LIG_IDX_LIST'][4] = 0
    data['LIG_GROUP'][4] = 0  # LIG_GROUP for 'pication

    data['FRAME'][5] = 5  # Add a new 'FRAME' value
    data['Prot_partner'][5] = 'A'  # Add a new 'Prot_partner' value
    data['INTERACTION'][5] = 'hbond'  # Add 'hbond' interaction
    data['ACCEPTORIDX'][5] = 301  # ACCEPTORIDX for 'hbond'
    data['DON_IDX'][5] = 0  # DON_IDX
    data['DONORTYPE'][5] = 0  # DON_IDX
    data['PROTISDON'][5] = True  # PROTISDON is True for 'hbond'
    data['ACCEPTOR_IDX'][5] = 0
    data['LIG_IDX_LIST'][5] = 0
    data['DONOR_IDX'][5] = 0
    data['LIG_GROUP'][5] = 0  # LIG_GROUP for 'pication

    # Add 'waterbridge' cases where PROTISDON is both True and False
    data['FRAME'][6] = 6  # Add a new 'FRAME' value
    data['Prot_partner'][6] = 'A'  # Add a new 'Prot_partner' value
    data['INTERACTION'][6] = 'waterbridge'  # Add 'waterbridge' interaction
    data['ACCEPTOR_IDX'][6] = 401  # ACCEPTOR_IDX for 'waterbridge'
    data['DON_IDX'][6] = 0  # DON_IDX
    data['DONORTYPE'][6] = 0  # DON_IDX
    data['DONOR_IDX'][6] = 0
    data['LIG_IDX_LIST'][6] = 0
    data['PROTISDON'][6] = True  # PROTISDON is True for 'waterbridge'
    data['LIG_GROUP'][6] = 0  # LIG_GROUP for 'pication

    data['FRAME'][7] = 7  # Add a new 'FRAME' value
    data['Prot_partner'][7] = 'B'  # Add a new 'Prot_partner' value
    data['INTERACTION'][7] = 'waterbridge'  # Add 'waterbridge' interaction
    data['DONOR_IDX'][7] = 501  # DONOR_IDX for 'waterbridge'
    data['DON_IDX'][7] = 0  # DON_IDX
    data['DONORTYPE'][7] = 0  # DON_IDX
    data['PROTISDON'][7] = False  # PROTISDON is False for 'waterbridge'
    data['ACCEPTOR_IDX'][7] = 0
    data['LIG_IDX_LIST'][7] = 0 # LIG_IDX_LIST for 'pication'
    data['LIG_GROUP'][7] = 0  # LIG_GROUP for 'pication

    # Add 'pistacking' case
    data['FRAME'][8] = 8  # Add a new 'FRAME' value
    data['Prot_partner'][8] = 'A'  # Add a new 'Prot_partner' value
    data['INTERACTION'][8] = 'pistacking'  # Add 'pistacking' interaction
    data['LIG_IDX_LIST'][8] = [7, 8]  # LIG_IDX_LIST for 'pistacking'
    data['LIG_GROUP'][8] = 0  # LIG_GROUP for 'pication
    data['ACCEPTOR_IDX'][8] = 0
    data['DON_IDX'][8] = 0  # DON_IDX
    data['DONOR_IDX'][8] = 0 
    data['PROTISDON'][8] = False
    data['DONORTYPE'][8] = 0  # DON_IDX

    # Add 'pication' case
    data['FRAME'][9] = 9  # Add a new 'FRAME' value
    data['Prot_partner'][9] = 'A'  # Add a new 'Prot_partner' value
    data['INTERACTION'][9] = 'pication'  # Add 'pication' interaction
    data['LIG_IDX_LIST'][9] = [9, 10]  # LIG_IDX_LIST for 'pication'
    data['LIG_GROUP'][9] = 'Group4'  # LIG_GROUP for 'pication'
    data['ACCEPTOR_IDX'][9] = 0
    data['DON_IDX'][9] = 0  # DON_IDX
    data['PROTISDON'][9] = False
    data['DONOR_IDX'][9] = 0 
    data['DONORTYPE'][9] = 0  # DON_IDX
    
    # Add 'metal' interaction case
    data['FRAME'][10] = 10  # Add a new 'FRAME' value
    data['Prot_partner'][10] = 'A'  # Add a new 'Prot_partner' value
    data['INTERACTION'][10] = 'metal'  # Add 'metal' interaction
    data['METAL_IDX'][10] = 401  # METAL_IDX for 'metal'
    data['METAL_TYPE'][10] = 'Fe'  # Metal type
    data['LOCATION'][10] = 'site1'  # Location
    data['ACCEPTOR_IDX'][10] = 0
    data['DONOR_IDX'][10] = 0

    data['FRAME'][11] = 11  # Add a new 'FRAME' value
    data['Prot_partner'][11] = 'A'  # Add a new 'Prot_partner' value
    data['INTERACTION'][11] = 'saltbridge'  # Add 'saltbridge' interaction
    data['LIG_IDX_LIST'][11] = [7, 8]  # Ligand index list for 'saltbridge PI'
    data['LIG_GROUP'][11] = 'Group4'  # Ligand group for 'saltbridge PI'
    data['PROTISPOS'][11] = False  # PROTISPOS is False for 'saltbridge PI'

    # Add 'hydrophobic' case where 'ring_found' is False
    data['FRAME'][12] = 12  # Add a new 'FRAME' value
    data['Prot_partner'][12] = 'C'  # Add a new 'Prot_partner' value
    data['INTERACTION'][12] = 'hydrophobic'  # Add 'hydrophobic' interaction
    data['LIGCARBONIDX'][12] = 104  # LIGCARBONIDX for 'hydrophobic' (not in any ring)

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
    3: {2: 'C_[5, 6]_Group3_NI_saltbridge'},
    4: {4: 'A_501_F_halogen'},
    5: {5: 'A_301_Acceptor_hbond'},
    6: {6: 'A_401_Acceptor_waterbridge'},
    7: {7: 'B_501_Donor_waterbridge'},
    8: {8: 'A_[7, 8]_pistacking'},
    9: {9: 'A_[9_ 10]_Group4_pication'},
    10: {10: 'A_401_Fe_site1_metal'},
    11: {11: 'A_[7, 8]_Group4_PI_saltbridge'},
    12: {12: 'C_104_hydrophobic'}
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

def test_unique_data_generation():
    # Test case 1: Check if the function returns an empty dictionary for an empty list
    result = unique_data_generation([])
    assert result == {}

    # Test case 2: Check if the function correctly identifies and stores unique values
    input_list = [1, 2, 2, 3, 3, 4, 5, 5]
    result = unique_data_generation(input_list)
    expected_result = {1: 1, 2: 2, 3: 3, 4: 4, 5: 5}
    assert result == expected_result

    # Test case 3: Check if the function handles strings
    input_list = ["apple", "banana", "apple", "cherry"]
    result = unique_data_generation(input_list)
    expected_result = {"apple": "apple", "banana": "banana", "cherry": "cherry"}
    assert result == expected_result

    
# Define a test case that uses the fixture
def test_remove_duplicate_values(test_remove_duplicates_data):
    input_data, expected_output = test_remove_duplicates_data
    assert remove_duplicate_values(input_data) == expected_output

def test_combine_subdict_values():
    # Test case 1: Empty input dictionary
    data = {}
    result = combine_subdict_values(data)
    assert result == {'all': []}

    # Test case 2: Input dictionary with sub-dictionaries
    data = {
        'dict1': {'a': 1, 'b': 2},
        'dict2': {'c': 3, 'd': 4},
        'dict3': {'e': 5, 'f': 6},
    }
    result = combine_subdict_values(data)
    assert result == {'all': [1, 2, 3, 4, 5, 6]}

    # Test case 3: Input dictionary with empty sub-dictionaries
    data = {
        'dict1': {},
        'dict2': {},
    }
    result = combine_subdict_values(data)
    assert result == {'all': []}

    # Test case 4: Input dictionary with sub-dictionaries containing various data types
    data = {
        'dict1': {'a': 1, 'b': 'text', 'c': [1, 2, 3]},
        'dict2': {'d': None, 'e': 5.5},
    }
    result = combine_subdict_values(data)
    assert result == {'all': [1, 'text', [1, 2, 3], None, 5.5]}

# Define a sample DataFrame for testing
sample_data = {
    'A': [1, 2, 3, 4, 5],
    'B': [2, 3, 4, 5, 6],
    'C': [3, 4, 5, 6, 7],
    'D': [4, 5, 6, 7, 8],
}
sample_df = pd.DataFrame(sample_data)

# Define the provided 'unique_columns_rings_grouped' data for testing
unique_columns_rings_grouped = {
    1: {0: 'A_101_hydrophobic'},
    2: {1: 'B_202_Donor_hbond', 3: 'A_102_hydrophobic'},
    3: {2: 'C_[5, 6]_Group3_NI_saltbridge'},
    4: {4: 'A_501_F_halogen'},
    5: {5: 'A_301_Acceptor_hbond'},
    6: {6: 'A_401_Acceptor_waterbridge'},
    7: {7: 'B_501_Donor_waterbridge'},
    8: {8: 'A_[7, 8]_pistacking'},
    9: {9: 'A_[9_ 10]_Group4_pication'},
    10: {10: 'A_401_Fe_site1_metal'},
    11: {11: 'A_[7, 8]_Group4_PI_saltbridge'},
    12: {12: 'C_104_hydrophobic'}
}

def test_filtering_values_with_provided_data():
    # Test case 1: Check if the function returns a list of filtered values
    threshold = 0.2  # 20% threshold
    frames = 1000  # Some arbitrary number of frames
    df = pd.DataFrame()  # Create an empty DataFrame for testing
    result = filtering_values(threshold, frames, df, unique_columns_rings_grouped)

    assert isinstance(result, list)

# Define a test case for the df_iteration_numbering function
def test_df_iteration_numbering():
    # Sample DataFrame for testing
    data = {
        'INTERACTION': ["hydrophobic", "hbond", "halogen"],
        'LIGCARBONIDX': [1, 2, 3],
        'Prot_partner': ['A', 'B', 'C'],
        'ACCEPTORIDX': [4, 5, 6],
        'PROTISDON': [True, False, True],
        'DONORIDX': [7, 8, 9],
        'DON_IDX': [10, 11, 12],
        'DONORTYPE': ['X', 'Y', 'Z'],
        'LIG_IDX_LIST': ['1,2', '3,4', '5,6'],
        'LIG_GROUP': ['Group1', 'Group2', 'Group3'],
    }
    df = pd.DataFrame(data)

    # Sample unique_data dictionary
    unique_data = {
        1: ["A_hydrophobic", "B_hydrophobic"],
        2: ["B_hbond_1_interaction", "C_hbond_2_interaction"],
        3: ["C_halogen_X_interaction", "D_halogen_Y_interaction"],
    }

    # Call the function with the sample DataFrame and unique_data
    df_iteration_numbering(df, unique_data)

    # Check if the DataFrame has been updated correctly
    assert df.loc[0, 'A_hydrophobic'] == 1
    assert df.loc[0, 'B_hydrophobic'] == 1
    assert df.loc[1, 'B_hbond_1_interaction'] == 1
    assert df.loc[1, 'C_hbond_2_interaction'] == 0
    assert df.loc[2, 'C_halogen_X_interaction'] == 0
    assert df.loc[2, 'D_halogen_Y_interaction'] == 0

