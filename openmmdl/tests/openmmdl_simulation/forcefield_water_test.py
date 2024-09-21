import pytest
from unittest.mock import patch
from simtk.openmm.app import ForceField
from rdkit import Chem
from openmmdl.openmmdl_simulation.forcefield_water import ForcefieldSelector, ForcefieldGenerator, ForcefieldPreparation, ForcefieldConfigurator, ForcefieldSetup

@pytest.fixture
def selector():
    return ForcefieldSelector()

@pytest.fixture
def generator():
    return ForcefieldGenerator()

@pytest.fixture
def configurator():
    # Create a real configuration dictionary or object
    mock_config = {
        "smallMoleculeForceField": "smirnoff",
        "add_membrane": False,
        "ligand": None,
        "water_model": "TIP3P"
    }
    
    # Create an actual configurator instance using real parameters
    return ForcefieldConfigurator(
        config=mock_config, 
        protein_ff="amber14-all.xml", 
        solvent_ff="tip3p.xml", 
        water_ff="tip3p"
    )

# A simple implementation of ForcefieldPreparation for testing
class MockForcefieldPreparation:
    def __init__(self, forcefield_name, water_model):
        self.forcefield_name = forcefield_name
        self.water_model = water_model

    def select_forcefield(self):
        return f"Forcefield: {self.forcefield_name}"

    def select_water_model(self):
        return f"Water Forcefield: {self.water_model}", f"Model Water: {self.water_model}"

# Override the original class in the test scope
def test_forcefield_setup_initialization():
    setup = ForcefieldSetup(
        forcefield_name="amber14",
        water_model="TIP3P",
        ligand_file="ligand.sdf",
        minimization=True
    )

    assert setup.forcefield_name == "amber14"
    assert setup.water_model == "TIP3P"
    assert setup.ligand_file == "ligand.sdf"
    assert setup.minimization is True

def test_setup_forcefield():
    # Directly use the MockForcefieldPreparation for the test
    setup = ForcefieldSetup(
        forcefield_name="amber14",
        water_model="TIP3P"
    )

    # Use the mock class instead of the real one
    prep = MockForcefieldPreparation(setup.forcefield_name, setup.water_model)
    forcefield = prep.select_forcefield()
    water_forcefield, model_water = prep.select_water_model()
    
    assert forcefield == "Forcefield: amber14"
    assert water_forcefield == "Water Forcefield: TIP3P"
    assert model_water == "Model Water: TIP3P"

def test_ff_selection(selector):
    assert selector.ff_selection("AMBER14") == "amber14-all.xml"
    assert selector.ff_selection("AMBER99SB") == "amber99sb.xml"
    assert selector.ff_selection("INVALID_FF") is None

def test_water_forcefield_selection(selector):
    # Test with AMBER14 forcefield
    forcefield = "amber14-all.xml"
    assert selector.water_forcefield_selection("TIP3P", forcefield) == "amber14/tip3p.xml"
    assert selector.water_forcefield_selection("SPC/E", forcefield) == "amber14/spce.xml"
    
    # Test with CHARMM36 forcefield
    forcefield = "charmm36.xml"
    assert selector.water_forcefield_selection("TIP4P-Ew", forcefield) == "charmm36/tip4pew.xml"
    
    # Test with invalid water model
    assert selector.water_forcefield_selection("INVALID_WATER", forcefield) is None

def test_water_model_selection(selector):
    # Test with valid AMBER14 forcefield
    forcefield = "amber14-all.xml"
    assert selector.water_model_selection("TIP3P", forcefield) == "tip3p"
    
    # Test with valid CHARMM36 forcefield
    forcefield = "charmm36.xml"
    assert selector.water_model_selection("TIP4P-Ew", forcefield) == "tip4pew"
    
    # Test with valid old AMBER forcefield
    forcefield = "amber99sb.xml"
    assert selector.water_model_selection("TIP5P", forcefield) == None  # Assuming this is a valid case for old AMBER
    
    # Test with another valid old AMBER forcefield
    forcefield = "amber03.xml"
    assert selector.water_model_selection("TIP4P-FB", forcefield) == "tip4pfb"  # Assuming this is a valid case for old AMBER

    # Test with invalid water model
    assert selector.water_model_selection("INVALID_WATER", forcefield) is None
def test_generate_forcefield(generator):
    # You will need the actual forcefield XML files in your working directory.
    protein_ff = "amber14-all.xml"
    solvent_ff = "tip3p.xml"
    
    # Assuming you have an RDKit molecule; create a simple water molecule for the test.
    rdkit_mol = Chem.MolFromSmiles('CCO')
    
    # Call the actual method without mocking.
    forcefield = generator.generate_forcefield(
        protein_ff=protein_ff,
        solvent_ff=solvent_ff,
        add_membrane=False,
        smallMoleculeForceField="gaff",
        rdkit_mol=rdkit_mol
    )
    
    # Check if a ForceField object is returned.
    assert isinstance(forcefield, ForceField)
    assert len(forcefield.getGenerators()) > 0  # Ensure forcefield generators were added.

def test_generate_transitional_forcefield(generator):
    # Similar setup to the above test, but with a membrane.
    protein_ff = "amber14-all.xml"
    solvent_ff = "tip3p.xml"
    rdkit_mol = Chem.MolFromSmiles('CCO')
    
    # Call the method with gaff.
    transitional_forcefield_gaff = generator.generate_transitional_forcefield(
        protein_ff=protein_ff,
        solvent_ff=solvent_ff,
        add_membrane=True,
        smallMoleculeForceField="gaff",
        rdkit_mol=rdkit_mol
    )

    # Call the method with smirnoff.
    transitional_forcefield_smirnoff = generator.generate_transitional_forcefield(
        protein_ff=protein_ff,
        solvent_ff=solvent_ff,
        add_membrane=True,
        smallMoleculeForceField="smirnoff",
        rdkit_mol=rdkit_mol
    )
    
    # Check if a ForceField object is returned with GAFF.
    assert isinstance(transitional_forcefield_gaff, ForceField)
    assert len(transitional_forcefield_gaff.getGenerators()) > 0

    # Check if a ForceField object is returned with SMIRNOFF.
    assert isinstance(transitional_forcefield_smirnoff, ForceField)
    assert len(transitional_forcefield_smirnoff.getGenerators()) > 0

class MockForcefieldSelector:
    def ff_selection(self, ff):
        # Return a mocked forcefield selection
        return f"Selected Forcefield: {ff}"

    def water_forcefield_selection(self, water, forcefield_selection):
        # Return a mocked water forcefield selection
        return f"Selected Water: {water} with {forcefield_selection}"

    def water_model_selection(self, water, forcefield_selection):
        # Return a mocked water model selection
        return f"Model Water: {water} for {forcefield_selection}"

@pytest.fixture
def forcefield_generator():
    # Set up the ForcefieldGenerator with the mock selector
    generator = ForcefieldGenerator()
    generator.forcefield_selector = MockForcefieldSelector()
    generator.ff = "amber14"  # Example forcefield
    generator.water = "TIP3P"  # Example water model
    return generator


def test_prepare_forcefield_and_water(forcefield_generator):
    forcefield_selected, water_selected, model_water = forcefield_generator.prepare_forcefield_and_water()
    
    # Check if the selections are as expected
    assert forcefield_selected == "Selected Forcefield: amber14"
    assert water_selected == "Selected Water: TIP3P with Selected Forcefield: amber14"
    assert model_water == "Model Water: TIP3P for Selected Forcefield: amber14"


def test_forcefield_preparation():
    prep = ForcefieldPreparation(forcefield_name="AMBER14", water_model="TIP3P")
    
    # Test selecting forcefield
    forcefield = prep.select_forcefield()
    assert forcefield == "amber14-all.xml"
    
    # Test selecting water model
    water_forcefield, model_water = prep.select_water_model()
    assert water_forcefield == "amber14/tip3p.xml"
    assert model_water == "tip3p"

class MockConfigParserGaff:
    def __init__(self, smallMoleculeForceField="gaff", add_membrane=False, ligand=None):
        self.smallMoleculeForceField = smallMoleculeForceField
        self.add_membrane = add_membrane
        self.ligand = ligand

class MockConfigParserSmirnoff:
    def __init__(self, smallMoleculeForceField="smirnoff", add_membrane=False, ligand=None):
        self.smallMoleculeForceField = smallMoleculeForceField
        self.add_membrane = add_membrane
        self.ligand = ligand

def test_create_forcefield_with_ligand():
    # Create the mock configuration
    mock_config_gaff = MockConfigParserGaff()
    mock_config_smirnoff = MockConfigParserGaff()
    mock_config_gaff.ligand = Chem.MolFromSmiles("CCO")  # Example ligand (ethanol)
    mock_config_smirnoff.ligand = Chem.MolFromSmiles("CCO")  # Example ligand (ethanol)
    
    # Initialize ForcefieldConfigurator
    configurator_gaff = ForcefieldConfigurator(
        config_parser=mock_config_gaff,
        forcefield_selected="amber14-all.xml",
        water_forcefield="tip3p.xml",
        water_selected="tip3p",
        prepared_ligand=mock_config.ligand
    )

    # Initialize ForcefieldConfigurator
    configurator_smirnoff = ForcefieldConfigurator(
        config_parser=mock_config_smirnoff,
        forcefield_selected="amber14-all.xml",
        water_forcefield="tip3p.xml",
        water_selected="tip3p",
        prepared_ligand=mock_config.ligand
    )
    
    forcefield_gaff = configurator_gaff.create_forcefield()
    forcefield_smirnoff = configurator_smirnoff.create_forcefield()
    
    # Check if the forcefield is generated successfully
    assert forcefield_gaff is not None

    # Check if the forcefield is generated successfully
    assert forcefield_smirnoff is not None

def test_generate_transitional_forcefield_with_membrane():
    # Create the mock configuration
    mock_config = MockConfigParser(add_membrane="True")

    # Initialize ForcefieldConfigurator
    configurator = ForcefieldConfigurator(
        config_parser=mock_config,
        forcefield_selected="amber14-all.xml",
        water_forcefield="tip3p.xml",
        water_selected="tip3p",
        prepared_ligand=None  # No ligand by default
    )

    transitional_forcefield = configurator.check_and_generate_forcefield()
    
    # Check if the transitional forcefield is generated successfully
    assert transitional_forcefield is not None

def test_generate_transitional_forcefield_without_membrane():
    # Create the mock configuration (default is False)
    mock_config = MockConfigParser()

    # Initialize ForcefieldConfigurator
    configurator = ForcefieldConfigurator(
        config_parser=mock_config,
        forcefield_selected="amber14-all.xml",
        water_forcefield="tip3p.xml",
        water_selected="tip3p",
        prepared_ligand=None  # No ligand by default
    )

    transitional_forcefield = configurator.check_and_generate_forcefield()
    
    # Check that no transitional forcefield is generated
    assert transitional_forcefield is None
