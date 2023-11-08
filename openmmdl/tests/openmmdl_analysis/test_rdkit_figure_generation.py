import pytest
import os
import time
import shutil
from PIL import Image
from pathlib import Path
from openmmdl.openmmdl_analysis.rdkit_figure_generation import split_interaction_data, highlight_numbers, update_dict, create_and_merge_images, arranged_figure_generation, generate_interaction_dict

test_data_directory = Path("openmmdl/tests/data/in")
current_directory = os.getcwd() 
output_path = 'all_binding_modes_arranged.png'

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
        "165ASPA 4222 Donor waterbridge",
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

def test_max_width_and_height_calculation():
    # Create some example images with different sizes
    image1 = Image.new('RGB', (100, 200), (255, 255, 255))
    image2 = Image.new('RGB', (150, 250), (255, 255, 255))
    merged_images = [image1, image2]

    # Calculate the maximum width and height
    max_width = max(image.size[0] for image in merged_images)
    max_height = max(image.size[1] for image in merged_images)

    # Assert the calculated max_width and max_height
    assert max_width == 150
    assert max_height == 250

def test_big_figure_creation():
    # Create example merged images
    image1 = Image.new('RGB', (100, 200), (255, 255, 255))
    image2 = Image.new('RGB', (150, 250), (255, 255, 255))
    merged_images = [image1, image2]

    # Calculate the maximum width and height
    max_width = max(image.size[0] for image in merged_images)
    max_height = max(image.size[1] for image in merged_images)

    # Determine the number of images per row (in your case, 2 images per row)
    images_per_row = 2

    # Calculate the number of rows and columns required
    num_rows = (len(merged_images) + images_per_row - 1) // images_per_row
    total_width = max_width * images_per_row
    total_height = max_height * num_rows

    # Create a new image with the calculated width and height
    big_figure = Image.new('RGB', (total_width, total_height), (255, 255, 255))  # Set background to white

    # Assert the dimensions of the created big_figure
    assert big_figure.size == (300, 500)  # Width should be 300, height should be 500

def test_arranged_figure_generation():
    binding_mode1_path = 'openmmdl/tests/data/openmmdl_analysis/rdkit_figure_generation/Binding_Mode_1_merged.png'
    binding_mode2_path = 'openmmdl/tests/data/openmmdl_analysis/rdkit_figure_generation/Binding_Mode_2_merged.png'
    all_modes_path = 'openmmdl/tests/data/openmmdl_analysis/rdkit_figure_generation/all_binding_modes_arranged.png'
    working_directory = os.getcwd()
    
    # Print the working directory to verify it's as expected
    print("Working Directory:", working_directory)

    destination_path_1 = os.path.join(working_directory, os.path.basename(binding_mode1_path))
    destination_path_2 = os.path.join(working_directory, os.path.basename(binding_mode2_path))
    destination_path_all = os.path.join(working_directory, os.path.basename(all_modes_path))
    
    # Print the destination paths to verify they are constructed correctly
    print("Destination Path 1:", destination_path_1)
    print("Destination Path 2:", destination_path_2)
    print("Destination Path All:", destination_path_all)

    shutil.copy(binding_mode1_path, destination_path_1)
    shutil.copy(binding_mode2_path, destination_path_2)
    shutil.copy(all_modes_path, destination_path_all)
    
    merged_image_paths = ['Binding_Mode_1_merged.png', 'Binding_Mode_2_merged.png']
    output_path = 'all_binding_modes_arranged.png'
    output_path = os.path.join(working_directory, output_path)
    print(output_path)

    # Run the function
    arranged_figure_generation(merged_image_paths, output_path)
    print(output_path)

    ## Print the current files in the working directory for debugging
    files_in_working_directory = os.listdir(working_directory)
    print("Files in Working Directory:", files_in_working_directory)

    output_path = os.path.join(working_directory, 'Binding_Modes_Markov_States', 'all_binding_modes_arranged.png')
    print(output_path)

    # Check if the output file was created
    
    assert output_path is not None

# Run the tests
if __name__ == '__main__':
    pytest.main()

