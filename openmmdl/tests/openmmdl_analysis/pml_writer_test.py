import numpy as np
import pandas as pd
import re
import os
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
import pytest
from openmmdl.openmmdl_analysis.pml_writer import *

#pml_writer tests
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

def test_generate_md_pharmacophore_cloudcenters(tmp_path):
    # Sample data for the DataFrame
    data = {
        'Acceptor_hbond_1': [1, 0, 1, 0, 1],
        'Donor_hbond_1': [0, 1, 0, 1, 0],
        'pistacking_1': [1, 0, 0, 1, 1],
        'hydrophobic_1': [0, 1, 0, 1, 0],
        'PI_saltbridge_1': [1, 0, 1, 0, 1],
        'NI_saltbridge_1': [0, 1, 0, 1, 0],
        'LIGCOO': ['(1.0, 2.0, 3.0)', '(2.0, 3.0, 4.0)', '(3.0, 4.0, 5.0)', '(4.0, 5.0, 6.0)', '(5.0, 6.0, 7.0)'],
        'PROTCOO': ['(7.0, 6.0, 5.0)', '(6.0, 5.0, 4.0)', '(5.0, 4.0, 3.0)', '(4.0, 3.0, 2.0)', '(3.0, 2.0, 1.0)'],
    }

    df = pd.DataFrame(data)

    # Output file paths
    output_filename = tmp_path / "test_output.pml"

    # Call the function
    generate_md_pharmacophore_cloudcenters(df, 'core_compound', output_filename, 'system_name', id_num=0)

    # Check if the output file is created
    assert os.path.isfile(output_filename), f"File {output_filename} not found."

    # Check if the generated XML is valid
    try:
        ET.parse(output_filename)
    except ET.ParseError:
        pytest.fail(f"Invalid XML in {output_filename}")

def test_generate_bindingmode_pharmacophore():
    dict_bindingmode = {
        "Acceptor_hbond_1": {
            "PROTCOO": [[7.0, 6.0, 5.0]],
            "LIGCOO": [[1.0, 2.0, 3.0]]
        },
        "hydrophobic_1": {
            "LIGCOO": [[4.0, 5.0, 6.0]]
        }
    }
    core_compound = "Ligand1"
    sysname = "System1"
    outname = "output_test"

    feature_types = {
        "Acceptor_hbond": "HBA",
        "Donor_hbond": "HBD",
        "pistacking": "AR",
        "hydrophobic": "H",
        "PI_saltbridge": "PI",
        "NI_saltbridge": "NI"
    }
    feature_id_counter = 0
    root = ET.Element("MolecularEnvironment", version="0.0", id=f"OpennMMDL_Analysis0", name=sysname)
    pharmacophore = ET.SubElement(root, "pharmacophore", name=sysname, id="pharmacophore0", pharmacophoreType="LIGAND_SCOUT")

    for interaction in dict_bindingmode.keys():
        # get feature type
        for interactiontype in feature_types.keys():
            if interactiontype in interaction:
                feature_type = feature_types[interactiontype]
                break
        # generate vector features
        if feature_type in ["HBA", "HBD"]:
            if feature_type == "HBA":
                orig_loc = dict_bindingmode[interaction]['PROTCOO'][0]
                targ_loc = dict_bindingmode[interaction]['LIGCOO'][0]
            elif feature_type == "HBD":
                orig_loc = dict_bindingmode[interaction]['LIGCOO'][0]
                targ_loc = dict_bindingmode[interaction]['PROTCOO'][0]
            feature_id_counter += 1
            points_to_lig = "true" if feature_type == "HBA" else "false"
            hasSyntheticProjectedPoint = "false"
            vector = ET.SubElement(
                pharmacophore,
                "vector",
                name=feature_type,
                featureId=interaction,
                pointsToLigand=points_to_lig,
                hasSyntheticProjectedPoint=hasSyntheticProjectedPoint,
                optional="false",
                disabled="false",
                weight="1.0",
                coreCompound=core_compound,
                id=f"feature{str(feature_id_counter)}"
            )
            origin = ET.SubElement(
                vector,
                "origin",
                x3=str(orig_loc[0]),
                y3=str(orig_loc[1]),
                z3=str(orig_loc[2]),
                tolerance="1.9499999"
            )
            target = ET.SubElement(
                vector,
                "target",
                x3=str(targ_loc[0]),
                y3=str(targ_loc[1]),
                z3=str(targ_loc[2]),
                tolerance="1.5"
            )
        # generate point features
        elif feature_type in ["H", "PI", "NI"]:
            position = dict_bindingmode[interaction]['LIGCOO'][0]
            feature_id_counter += 1
            point = ET.SubElement(
                pharmacophore,
                "point",
                name=feature_type,
                featureId=interaction,
                optional="false",
                disabled="false",
                weight="1.0",
                coreCompound=core_compound,
                id=f"feature{str(feature_id_counter)}"
            )
            position = ET.SubElement(
                point,
                "position",
                x3=str(position[0]),
                y3=str(position[1]),
                z3=str(position[2]),
                tolerance="1.5"
            )
        # generate plane features
        elif feature_type == "AR":
            feature_id_counter += 1
            lig_loc = dict_bindingmode[interaction]['LIGCOO'][0]
            prot_loc = dict_bindingmode[interaction]['PROTCOO'][0]

            # normalize vector of plane
            vector = np.array(lig_loc) - np.array(prot_loc)
            normal_vector = vector / np.linalg.norm(vector)
            x, y, z = normal_vector

            plane = ET.SubElement(pharmacophore,
                                  "plane",
                                  name=feature_type,
                                  featureId=interaction,
                                  optional="false",
                                  disabled="false",
                                  weight="1.0",
                                  coreCompound=core_compound,
                                  id=f"feature{str(feature_id_counter)}")
            position = ET.SubElement(plane,
                                     "position",
                                     x3=str(lig_loc[0]),
                                     y3=str(lig_loc[1]),
                                     z3=str(lig_loc[2]),
                                     tolerance="0.9")
            normal = ET.SubElement(plane,
                                   "normal",
                                   x3=str(x),
                                   y3=str(y),
                                   z3=str(z),
                                   tolerance="0.43633232")

    tree = ET.ElementTree(root)
    tree_str = ET.tostring(root, encoding="unicode")
    print(f"XML Structure:\n{tree_str}")

    hydrophobic_point = root.find(".//point[@name='hydrophobic']")
    print(f"Found hydrophobic_point: {hydrophobic_point}")
    assert hydrophobic_point is not None, "Hydrophobic point not found"



def test_generate_pharmacophore_centers_all_points():
    # Sample data for the DataFrame
    data = {
        'interaction1': [1, 0, 1, 0, 1],
        'interaction2': [0, 1, 0, 1, 0],
        'LIGCOO': ['(1.0, 2.0, 3.0)', '(2.0, 3.0, 4.0)', '(3.0, 4.0, 5.0)', '(4.0, 5.0, 6.0)', '(5.0, 6.0, 7.0)'],
    }

    df = pd.DataFrame(data)

    # Sample interactions
    interactions = ['interaction1', 'interaction2']

    # Call the function
    pharmacophore = generate_pharmacophore_centers_all_points(df, interactions)

    # Check if the generated pharmacophore has the expected structure
    assert isinstance(pharmacophore, dict), "Pharmacophore should be a dictionary."
    
    for interaction in interactions:
        assert interaction in pharmacophore, f"{interaction} not found in the generated pharmacophore."

        points = pharmacophore[interaction]
        assert isinstance(points, list), f"Pharmacophore points for {interaction} should be a list."
        
        # Check if the points have the expected structure
        for point in points:
            assert isinstance(point, list) and len(point) == 3, "Each point should be a list of three coordinates."
    


def test_generate_point_cloud_pml(tmp_path):
    # Sample data for the cloud_dict
    cloud_dict = {
        'feature1': {
            'interaction1': [(1.0, 2.0, 3.0), (1.5, 2.5, 3.5), (2.0, 3.0, 4.0)],
            'interaction2': [(2.0, 3.0, 4.0), (2.5, 3.5, 4.5), (3.0, 4.0, 5.0)],
        },
        'feature2': {
            'interaction3': [(3.0, 4.0, 5.0), (3.5, 4.5, 5.5), (4.0, 5.0, 6.0)],
        },
    }

    # Output file paths
    outname = tmp_path / "test_output"
    outname_pml = f"{outname}.pml"

    # Call the function
    generate_point_cloud_pml(cloud_dict, "system_name", outname)

    # Check if the output file is created
    assert os.path.isfile(outname_pml), f"File {outname_pml} not found."

    # Check if the generated XML is valid
    try:
        ET.parse(outname_pml)
    except ET.ParseError:
        pytest.fail(f"Invalid XML in {outname_pml}")
