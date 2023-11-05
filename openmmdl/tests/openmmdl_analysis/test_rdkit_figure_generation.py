import pytest
import os
import time
import shutil
from PIL import Image
from pathlib import Path
from openmmdl.openmmdl_analysis.rdkit_figure_generation import split_interaction_data, highlight_numbers, update_dict, create_and_merge_images, arranged_figure_generation

test_data_directory = Path("openmmdl/tests/data/in")
current_directory = os.getcwd() 

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


# Define test data
@pytest.fixture
def merged_image_paths():
    image_paths = []
    for i in range(1, 5):
        image = Image.new('RGB', (100, 100), (i * 25, i * 25, i * 25))
        image_path = f"{current_directory}/image_{i}.png"
        image.save(image_path)
        image_paths.append(image_path)
    return image_paths

@pytest.fixture
def output_path():
    return (f"{current_directory}/output.png")

# Test the arranged_figure_generation function
def test_arranged_figure_generation(merged_image_paths,output_path):
    arranged_figure_generation(merged_image_paths, output_path)
    assert os.path.exists(output_path)
    
    # Check if the output file is an image
    with Image.open(output_path) as output_image:
        assert output_image.mode == 'RGB'

    # Check the output image dimensions (you may need to adjust this depending on your input)
    with Image.open(output_path) as output_image:
        expected_width = 200  # 2 images per row
        expected_height = 200  # 2 rows
        assert output_image.size == (expected_width, expected_height)

    # Check if individual image files are removed
    for path in merged_image_paths:
        assert not os.path.exists(path)

    # Check if the output file is renamed
    new_output_path = "Binding_Modes_Markov_States/output.png"
    assert os.path.exists(new_output_path)



# Run the tests
if __name__ == '__main__':
    pytest.main()

