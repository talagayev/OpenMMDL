import pytest
import os
from pathlib import Path
from openmmdl.openmmdl_analysis.rdkit_figure_generation import split_interaction_data


@pytest.mark.parametrize("input_data, expected_output", [
    (["60GLUA_4206_4207_4216_4217_4218_4205_hydrophobic"], ['60GLUA 4206 4207 4216 4217 4218 4205 hydrophobic']),
    (["165ASPA_4203_Acceptor_hbond"], ['165ASPA 4203 Acceptor hbond']),
    (["125TYRA_4192_Acceptor_waterbridge"], ['125TYRA 4192 Acceptor waterbridge']),
])
def test_split_interaction_data(input_data, expected_output):
    result = split_interaction_data(input_data)
    assert result == expected_output


# Define test cases
@pytest.mark.parametrize("split_data, starting_idx, expected_output_highlight", [
    (["Res1_TypeA 1 2 3 HydrogenBond", "Res2_TypeB 4 5 6 HydrogenBond"], 1, ([0], [1, 2], [], [], [], [], [], [], [], [], [])),
    (["Res1_TypeA 1 2 3 HydrogenBond Donor", "Res2_TypeB 4 5 6 HydrogenBond Acceptor"], 1, ([0], [1], [2], [], [], [], [], [], [], [], [])),
    (["Res1_TypeA 1 2 3 HydrogenBond Donor", "Res2_TypeB 4 5 6 HydrogenBond Donor"], 1, ([0, 1], [], [], [], [], [], [], [], [], [], [])),
    (["Res1_TypeA 1 2 Hydrophobic", "Res2_TypeB 4 5 Hydrophobic"], 1, ([], [], [], [0, 1], [], [], [], [], [], [], [])),
    (["Res1_TypeA 1 2 3 WaterBridge", "Res2_TypeB 4 5 WaterBridge"], 1, ([], [], [], [], [0, 1], [], [], [], [], [], [])),
    (["Res1_TypeA 1 2 3 Pistacking", "Res2_TypeB 4 5 6 Pistacking"], 1, ([], [], [], [], [], [0, 1, 2], [], [], [], [], [])),
    (["Res1_TypeA 1 2 3 Halogen HalogenType", "Res2_TypeB 4 5 6 Halogen HalogenType"], 1, ([], [], [], [], [], [], [0, 1, 2], [], [], [], [])),
    (["Res1_TypeA 1,2 3 SaltBridge PI", "Res2_TypeB 4,5 SaltBridge NI"], 1, ([], [], [], [], [], [], [], [0, 1], [2, 3], [], [])),
    (["Res1_TypeA 1 2 Pication HalogenType", "Res2_TypeB 4 5 Pication HalogenType"], 1, ([], [], [], [], [], [], [], [], [], [0, 1, 2], [])),
    (["Res1_TypeA 1 2 3 Metal HalogenType", "Res2_TypeB 4 5 6 Metal HalogenType"], 1, ([], [], [], [], [], [], [], [], [], [], [0, 1, 2, 3])),
])
def test_highlight_numbers(split_data, starting_idx, expected_output_highlight):
    result = highlight_numbers(split_data, starting_idx)
    assert result == expected_output_highlight

# Run the tests
if __name__ == "__main__":
    pytest.main()
