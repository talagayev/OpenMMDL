import pytest
import os
import time
import shutil
from PIL import Image
from pathlib import Path
from openmmdl.openmmdl_analysis.rdkit_figure_generation import split_interaction_data, highlight_numbers, update_dict, create_and_merge_images, arranged_figure_generation, generate_interaction_dict

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

def test_generate_interaction_dict():
    # Test with a known interaction type 'hydrophobic'
    interaction_type = 'hydrophobic'
    keys = [1, 2, 3]
    expected_result = {
        1: (1.0, 1.0, 0.0),
        2: (1.0, 1.0, 0.0),
        3: (1.0, 1.0, 0.0)
    }
    result = generate_interaction_dict(interaction_type, keys)
    assert result == expected_result

# Define test data
binding_mode = "Binding_Mode_1"
occurrence_percent = 92.0
split_data = [
    '161PHEA 4221 Acceptor hbond',
    'FRAME  FRAME',
    '207ILEA 4205 4206 4207 4208 4209 4204 hydrophobic',
    '166ARGA 4220,4221 Carboxylate NI saltbridge'
]
merged_image_paths = []

# Create a fixture to prepare any necessary resources (e.g., test images)
@pytest.fixture
def prepare_resources():
    # Copy the image from the specified path to the working directory
    working_directory = os.getcwd()
    existing_image_path = "openmmdl/tests/data/in/Binding_Mode_1.png"
    copied_image_path = os.path.join(working_directory, "Binding_Mode_1.png")
    shutil.copy(existing_image_path, copied_image_path)

# Define the test function
def test_create_and_merge_images(prepare_resources):
    # Path to the existing image in the working directory
    existing_image_path = "Binding_Mode_1.png"

    # Call the function with the test data
    merged_image_paths = create_and_merge_images(binding_mode, occurrence_percent, split_data)

    # Assert that the function returns a list of image paths
    assert isinstance(merged_image_paths, list)
    assert all(isinstance(path, str) for path in merged_image_paths)

    # Assert that the created images exist and can be opened
    for image_path in merged_image_paths:
        with Image.open(image_path) as img:
            assert img is not None

    # Assert that the existing image can be opened
    with Image.open(existing_image_path) as existing_img:
        assert existing_img is not None


# Run the tests
if __name__ == '__main__':
    pytest.main()

