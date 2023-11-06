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
        'PROTISPOS': {0: True, 1: False, 2: True, 3: True},
        'DON_IDX': {0: 0, 1: 0, 2: 0, 3: 0},
        'DONORTYPE': {0: 0, 1: 0, 2: 0, 3: 0},
        'ACCEPTOR_IDX': {0: 0, 1: 0, 2: 0, 3: 0},
        'DONOR_IDX': {0: 0, 1: 0, 2: 0, 3: 0},
        'LIG_GROUP': {0: 0, 1: 0, 2: 0, 3: 0},
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
    9: {9: 'A_[9, 10]_Group4_pication'}
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
    
