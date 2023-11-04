import pytest
import os
import shutil
from pathlib import Path
from openmmdl.openmmdl_analysis.rdkit_figure_generation import split_interaction_data, highlight_numbers, update_dict, create_and_merge_images

test_data_directory = Path("openmmdl/tests/data/in")

@pytest.mark.parametrize("input_data, expected_output", [
    (["60GLUA_4206_4207_4216_4217_4218_4205_hydrophobic"], ['60GLUA 4206 4207 4216 4217 4218 4205 hydrophobic']),
    (["165ASPA_4203_Acceptor_hbond"], ['165ASPA 4203 Acceptor hbond']),
    (["125TYRA_4192_Acceptor_waterbridge"], ['125TYRA 4192 Acceptor waterbridge']),
])
def test_split_interaction_data(input_data, expected_output):
    result = split_interaction_data(input_data)
    assert result == expected_output

def test_highlight_numbers():
    # Input data
    split_data = [
        "163GLYA 4202 Acceptor hbond",
        "165ASPA 4203 Donor hbond",
        "165ASPA 4222 Donor hbond",
        "165ASPA 4203 Acceptor hbond",
        "125TYRA 4192 Acceptor waterbridge",
        "161PHEA 4211 4212 4213 4214 4215 4210 hydrophobic",
        "59ARGA 4205 4206 4207 4216 4217 4218 Aromatic pication",
        "59ARGA 4194 F halogen",
        "166ARGA 4202,4203 Carboxylate NI saltbridge"
    ]

    starting_idx = 1  # Updated starting index

    result = highlight_numbers(split_data, starting_idx)

    highlighted_hbond_donor, highlighted_hbond_acceptor, highlighted_hbond_both, \
    highlighted_hydrophobic, highlighted_waterbridge, highlighted_pistacking, highlighted_halogen, \
    highlighted_ni, highlighted_pi, highlighted_pication, highlighted_metal = result

    assert highlighted_hbond_donor is not None
    assert highlighted_hbond_acceptor is not None
    assert highlighted_hbond_both is not None
    assert highlighted_hydrophobic is not None
    assert highlighted_waterbridge is not None
    assert highlighted_halogen is not None
    assert highlighted_ni is not None
    assert highlighted_pication is not None
    
def test_update_dict():
    # Test case 1: Check if the target dictionary is updated correctly
    target_dict = {1: '1', 2: '2'}
    source_dict = {3: '3', 4: '4'}
    update_dict(target_dict, source_dict)
    assert target_dict == {1: '1', 2: '2', 3: '3', 4: '4'}

    # Test case 2: Check if the function handles multiple source dictionaries
    target_dict = {}
    source_dict1 = {1: '1'}
    source_dict2 = {2: '2', 3: '3'}
    update_dict(target_dict, source_dict1, source_dict2)
    assert target_dict == {1: '1', 2: '2', 3: '3'}

    # Test case 3: Check if the function handles empty source dictionaries
    target_dict = {1: '1', 2: '2'}
    update_dict(target_dict)  # No source dictionaries provided
    assert target_dict == {1: '1', 2: '2'}

# Define test cases
@pytest.mark.parametrize("binding_mode, occurrence_percent, split_data, expected_output", [
    ("Binding_mode_1", {"Binding_mode_1": 60}, ["163GLYA 4202 Acceptor hbond", "165ASPA 4222 Donor hbond", "161PHEA 4211 4212 4213 4214 4215 4210 hydrophobic"], ["Binding_mode_1_merged.png"]),
    ("Binding_mode_2", {"Binding_mode_2": 20}, ["59ARGA 4194 F halogen", "125TYRA 4192 Acceptor waterbridge", "166ARGA 4202,4203 Carboxylate NI saltbridge"], ["Binding_mode_2_merged.png"]),
])
def test_create_and_merge_images(binding_mode, occurrence_percent, split_data, expected_output):

    original_cwd = os.getcwd()
    binding_mode_1 = f"{test_data_directory}/Binding_Mode_1.png"
    binding_mode_2 = f"{test_data_directory}/Binding_Mode_2.png"
    binding_mode_1 = os.path.join(test_data_directory, "Binding_Mode_1.png")
    binding_mode_2 = os.path.join(test_data_directory, "Binding_Mode_2.png")
    shutil.copy(binding_mode_1, original_cwd)
    shutil.copy(binding_mode_2, original_cwd)
        
    merged_image_paths = create_and_merge_images(binding_mode, occurrence_percent, split_data, [])
    
    assert "merged.png" in merged_image_paths
    for merged_image_path in merged_image_paths:
        assert os.path.exists(merged_image_path)

# Run the tests
if __name__ == '__main__':
    pytest.main()



# Run the tests
if __name__ == "__main__":
    pytest.main()
