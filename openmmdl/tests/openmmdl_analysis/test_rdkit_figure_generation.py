import pytest
import os
from pathlib import Path
from openmmdl.openmmdl_analysis.rdkit_figure_generation import split_interaction_data, highlight_numbers


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
    

# Run the tests
if __name__ == "__main__":
    pytest.main()
