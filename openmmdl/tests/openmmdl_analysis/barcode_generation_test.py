import numpy as np
import pandas as pd
import re
import os
import matplotlib.pyplot as plt
import pytest
from openmmdl.openmmdl_analysis.barcode_generation import *

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

def test_plot_barcodes(tmp_path):
    # Test case 1: No barcodes provided
    with pytest.raises(SystemExit):
        plot_barcodes([], "no_barcodes.png")

    # Test case 2: Single barcode
    barcode = np.random.randint(0, 2, 100)
    save_path = tmp_path / "single_barcode.png"
    plot_barcodes([("Barcode 1", barcode)], save_path)

    assert save_path.is_file()
    
    # Test case 3: Multiple barcodes
    barcodes = [("Barcode 1", np.random.randint(0, 2, 100)),
                ("Barcode 2", np.random.randint(0, 2, 100)),
                ("Barcode 3", np.random.randint(0, 2, 100))]
    save_path = tmp_path / "multiple_barcodes.png"
    plot_barcodes(barcodes, save_path)
    
    assert save_path.is_file()

