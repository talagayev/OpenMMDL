import os
import pytest
from pathlib import Path
import pandas as pd
import numpy as np
import mdtraj as md
from plip.structure.preparation import PDBComplex, LigandFinder, Mol, PLInteraction, create_df_from_binding_site

from openmmdl.openmmdl_analysis.interaction_gathering import characterize_complex, retrieve_plip_interactions


test_data_directory = Path("openmmdl/tests/data/in")
topology_file = f"{test_data_directory}/complex.pdb"
binding_site_id = "UNK:X:0"
lig_name = "UNK"

# Test the function
def test_characterize_complex():
    # Call the function
    interaction_set = characterize_complex(topology_file, binding_site_id)

    # Check if the function returns a PLInteraction object
    assert isinstance(interaction_set, PLInteraction)

def test_retrieve_plip_interactions():
    # Call the function
    interactions = retrieve_plip_interactions(topology_file, lig_name)

    # Check if the function returns a dictionary
    assert isinstance(interactions, dict)

# Define test data
sample_interactions = {
    "hydrophobic": [["Column1", "Column2"], [1, 2], [3, 4]],
    "hbond": [["ColumnA", "ColumnB"], ["A", "B"], ["C", "D"]],
}

def test_create_df_from_binding_site():
    # Test with valid interaction type
    df = create_df_from_binding_site(sample_interactions, interaction_type="hydrophobic")
    assert isinstance(df, pd.DataFrame)
    assert df.shape == (2, 2)
    assert list(df.columns) == ["Column1", "Column2"]

    # Test with default interaction type
    df_default = create_df_from_binding_site(sample_interactions)
    assert isinstance(df_default, pd.DataFrame)
    assert df_default.shape == (2, 2)
    assert list(df_default.columns) == ["ColumnA", "ColumnB"]

    # Test with an invalid interaction type (should default to 'hbond')
    df_invalid = create_df_from_binding_site(sample_interactions, interaction_type="invalid_type")
    assert isinstance(df_invalid, pd.DataFrame)
    assert df_invalid.shape == (2, 2)
    assert list(df_invalid.columns) == ["ColumnA", "ColumnB"]

if __name__ == "__main":
    pytest.main()
