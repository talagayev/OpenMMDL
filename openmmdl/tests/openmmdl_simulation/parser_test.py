import pytest
import os
from openmmdl.openmmdl_simulation.parser import ConfigParser  # Adjust the import based on your file structure

@pytest.fixture
def temp_config_file(tmp_path):
    # Create a temporary config file for testing
    config_data = """
    VAR1 = 'value1'
    VAR2 = "value2"
    VAR3 = value3
    """
    config_file = tmp_path / "config.txt"
    config_file.write_text(config_data)
    return config_file

def test_initialization_and_parsing(temp_config_file):
    parser = ConfigParser(temp_config_file)
    assert parser.variables == {
        'VAR1': 'value1',
        'VAR2': 'value2',
        'VAR3': 'value3'
    }

def test_get_variable_existing(temp_config_file):
    parser = ConfigParser(temp_config_file)
    assert parser.get_variable('VAR1') == 'value1'

def test_get_variable_non_existing(temp_config_file):
    parser = ConfigParser(temp_config_file)
    assert parser.get_variable('NON_EXISTENT', default='default_value') == 'default_value'

def test_attribute_access(temp_config_file):
    parser = ConfigParser(temp_config_file)
    assert parser.VAR1 == 'value1'
    assert parser.VAR2 == 'value2'
    assert parser.VAR3 == 'value3'
    assert parser.NON_EXISTENT is None  # Should not raise an error

def test_print_variables(capsys, temp_config_file):
    parser = ConfigParser(temp_config_file)
    parser.print_variables()
    captured = capsys.readouterr()
    expected_output = "VAR1 = value1\nVAR2 = value2\nVAR3 = value3\n"
    assert captured.out == expected_output
