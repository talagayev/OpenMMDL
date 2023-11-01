import pytest
import simtk.openmm.app as app
from openff.toolkit.topology import Molecule
from openmmforcefields.generators import GAFFTemplateGenerator
from openmmdl.openmmdl_simulation.scripts.forcefield_water import ff_selection, water_forcefield_selection, water_model_selection, generate_forcefield, generate_transitional_forcefield

# Replace 'your_module' with the actual name of the module containing your functions.

@pytest.fixture
def sample_rdkit_molecule():
    """
    A sample RDKit molecule for testing.
    """
    from rdkit import Chem
    mol = Chem.MolFromSmiles("CCO")
    return mol

def test_ff_selection():
    assert ff_selection('AMBER14') == 'amber14-all.xml'
    assert ff_selection('CHARMM36') == 'charmm36.xml'
    assert ff_selection('NonexistentFF') is None

def test_water_forcefield_selection():
    assert water_forcefield_selection('TIP3P', 'amber14-all.xml') == 'amber14/tip3p.xml'
    assert water_forcefield_selection('SPC/E', 'charmm36.xml') == 'charmm36/spce.xml'
    assert water_forcefield_selection('TIP3P', 'NonexistentFF') is None

def test_water_model_selection():
    assert water_model_selection('TIP3P', 'amber14-all.xml') == 'tip3p'
    assert water_model_selection('TIP3P-PME-B', 'charmm36.xml') == 'charmm'
    assert water_model_selection('TIP3P', 'NonexistentFF') is None

def test_generate_forcefield(sample_rdkit_molecule):
    forcefield = generate_forcefield('amber14-all.xml', 'amber14/tip3p.xml', False, sample_rdkit_molecule)
    assert isinstance(forcefield, app.ForceField)

@pytest.mark.parametrize("add_membrane, protein_ff, expected_ff", [
    (True, 'amber10.xml', ['amber10.xml', 'amber14/lipid17.xml']),
    (True, 'amber99sb.xml', ['amber99sb.xml', 'amber14/lipid17.xml']),
    (False, 'amber14-all.xml', ['amber14-all.xml', 'amber14/tip3p.xml']),
    (False, 'amber03.xml', ['amber03.xml', 'amber14/tip3p.xml']),
])
def test_generate_forcefield_membrane_logic(add_membrane, protein_ff, expected_ff):
    forcefield = generate_forcefield(protein_ff, 'amber14/tip3p.xml', add_membrane, sample_rdkit_molecule)
    ff_files = [ff._forcefield for ff in forcefield._forces]
    assert ff_files == expected_ff

def test_generate_transitional_forcefield(sample_rdkit_molecule):
    transitional_forcefield = generate_transitional_forcefield('amber14-all.xml', 'tip3p.xml', True, sample_rdkit_molecule)
    assert isinstance(transitional_forcefield, app.ForceField)
