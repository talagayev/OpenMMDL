import pytest
from openmmdl.openmmdl_analysis.ligand_processing import increase_ring_indices, convert_ligand_to_smiles


def test_increase_ring_indices():
    # Test case 1: Check if ring indices are correctly increased
    ring = [1, 2, 3]
    lig_index = 10
    result = increase_ring_indices(ring, lig_index)
    assert result == [11, 12, 13]

    # Test case 2: Check with a different lig_index
    ring = [3, 4, 5]
    lig_index = 20
    result = increase_ring_indices(ring, lig_index)
    assert result == [23, 24, 25]

    # Add more test cases as needed

if __name__ == '__main__':
    pytest.main()
