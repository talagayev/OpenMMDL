import pytest
import openmmdl
from openmmdl.openmmdl_setup.openmmdlsetup import configureDefaultAmberOptions

@pytest.fixture
def session():
    # Implement a fixture that provides a session-like dictionary
    return {}

def test_configure_default_amber_options():
    # Call the function without any arguments
    configureDefaultAmberOptions()

    # Check that the default options are correctly set
    # Replace the following checks with your specific expectations
    assert session['nmLig'] is False
    assert session['lig_ff'] == 'gaff'
    assert session['charge_value'] == '0'
    assert session['charge_method'] == 'bcc'
    assert session['spLig'] is False
    assert session['prot_ff'] == 'ff14SB'
    assert session['dna_ff'] == 'OL15'
    assert session['rna_ff'] == 'OL3'
    assert session['carbo_ff'] == 'GLYCAM_06j'
    assert session['addType'] == 'addWater'
    assert session['boxType'] == 'cube'
    assert session['dist'] == '10'
    assert session['lipid_tp'] == 'POPC'
    assert session['other_lipid_tp_input'] == 'POPC:TOPC'
    assert session['lipid_ratio'] == '1:1'
    assert session['lipid_ff'] == 'lipid21'
    assert session['dist2Border'] == '15'
    assert session['padDist'] == '17'
    assert session['water_ff'] == 'tip3p'
    assert session['pos_ion'] == 'Na+'
    assert session['neg_ion'] == 'Cl-'
    assert session['ionConc'] == '0.15'
