import pytest
import os
from pathlib import Path
from openmmdl.openmmdl_analysis.rdkit_figure_generation import split_interaction_data


@pytest.mark.parametrize("input_data, expected_output", [
    (["Res1_TypeA_1 2 3_HydrogenBond"], ['Res1_TypeA 1 2 3 HydrogenBond']),
    (["Res2_TypeB_4 5 6_Electrostatic"], ['Res2_TypeB 4 5 6 Electrostatic']),
    (["Res3_TypeC_7 8_InteractionType"], ['Res3_TypeC 7 8 InteractionType']),
    (["Res4_TypeD_9 10_Unknown"], ['Res4_TypeD 9 10 Unknown']),
    (["Res5_TypeE_11 12 13_AnotherType"], ['Res5_TypeE 11 12 13 AnotherType']),
])
def test_split_interaction_data(input_data, expected_output):
    result = split_interaction_data(input_data)
    assert result == expected_output

# Run the tests
if __name__ == "__main__":
    pytest.main()
