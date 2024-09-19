import pytest
from openmmdl.openmmdl_setup.amberscript_creator import AmberScriptGenerator  # replace 'your_module' with the actual module name

@pytest.fixture
def mock_data():
    """Fixture to provide mock session and uploadedFiles data."""
    return {
        "session": {
            "rcpType": "protRcp",
            "prot_ff": "ff14SB",
            "other_prot_ff_input": "other_prot_ff_input",
            "dna_ff": "ff99bsc1",
            "other_dna_ff_input": "other_dna_ff_input",
            "rna_ff": "ff99",
            "other_rna_ff_input": "other_rna_ff_input",
            "carbo_ff": "glycam",
            "other_carbo_ff_input": "other_carbo_ff_input",
            "charge_method": "bcc",
            "charge_value": "-1",
            "lig_ff": "gaff",
            "nmLig": True,
            "spLig": True,
            "addType": "addWater",
            "boxType": "cube",
            "dist": "10.0"
        },
        "uploadedFiles": {
            "protFile": [("file1", "protein.pdb")],
            "dnaFile": [("file2", "dna.pdb")],
            "rnaFile": [("file3", "rna.pdb")],
            "carboFile": [("file4", "carbo.pdb")],
            "nmLigFile": [("file5", "ligand.pdb")],
            "spLigFile": [("file6", "ligand.pdb")],
            "prepcFile": [("file7", "ligand.prepc")],
            "frcmodFile": [("file8", "ligand.frcmod")]
        }
    }


@pytest.fixture
def base_mock_data():
    """Fixture providing mock data for different receptor types."""
    return {
        "protRcp": {
            "session": {
                "rcpType": "protRcp",
                "prot_ff": "ff14SB",
                "other_prot_ff_input": "custom_ff"
            },
            "uploadedFiles": {
                "protFile": [["file1", "protein.pdb"]]
            }
        },
        "dnaRcp": {
            "session": {
                "rcpType": "dnaRcp",
                "dna_ff": "bsc1",
                "other_dna_ff_input": "custom_dna_ff"
            },
            "uploadedFiles": {
                "dnaFile": [["file2", "dna.pdb"]]
            }
        },
        "rnaRcp": {
            "session": {
                "rcpType": "rnaRcp",
                "rna_ff": "ff99SB",
                "other_rna_ff_input": "custom_rna_ff"
            },
            "uploadedFiles": {
                "rnaFile": [["file3", "rna.pdb"]]
            }
        },
        "carboRcp": {
            "session": {
                "rcpType": "carboRcp",
                "carbo_ff": "GLYCAM",
                "other_carbo_ff_input": "custom_carbo_ff"
            },
            "uploadedFiles": {
                "carboFile": [["file4", "carbo.pdb"]]
            }
        }
    }

def test_add_openmmdl_logo(mock_data):
    """Test if add_openmmdl_logo correctly appends the logo to the amber_script list."""
    session, uploadedFiles = mock_data
    amber_script_gen = AmberScriptGenerator(session, uploadedFiles)
    
    # Prepare an empty amber_script list
    amber_script = []

    # Call the method
    amber_script_gen.add_openmmdl_logo(amber_script)

    # Define the expected logo output
    expected_logo = """
#       ,-----.    .-------.     .-''-.  ,---.   .--.,---.    ,---.,---.    ,---. ______       .---.      
#     .'  .-,  '.  \\  _(`)_ \\  .'_ _   \\ |    \\  |  ||    \\  /    ||    \\  /    ||    _ `''.   | ,_|      
#    / ,-.|  \\ _ \\ | (_ o._)| / ( ` )   '|  ,  \\ |  ||  ,  \\/  ,  ||  ,  \\/  ,  || _ | ) _  \\,-./  )      
#   ;  \\  '_ /  | :|  (_,_) /. (_ o _)  ||  |\\_ \\|  ||  |\\_   /|  ||  |\\_   /|  ||( ''_'  ) |\\  '_ '`)    
#   |  _`,/ \\ _/  ||   '-.-' |  (_,_)___||  _( )_\\  ||  _( )_/ |  ||  _( )_/ |  || . (_) `. | > (_)  )    
#   : (  '\\_/ \\   ;|   |     '  \\   .---.| (_ o _)  || (_ o _) |  || (_ o _) |  ||(_    ._) '(  .  .-'    
#    \\ `"/  \\  ) / |   |      \\  `-'    /|  (_,_)\  ||  (_,_)  |  ||  (_,_)  |  ||  (_.\\.' /  `-'`-'|___  
#     '. \\_/``".'  /   )       \\       / |  |    |  ||  |      |  ||  |      |  ||       .'    |        \\ 
#       '-----'    `---'        `'-..-'  '--'    '--''--'      '--''--'      '--''-----'`      `--------`                                                 
        """

    # Verify the logo is appended
    assert len(amber_script) == 1
    assert amber_script[0] == expected_logo

def test_add_openmmdl_amber_logo(mock_data):
    """Test if add_openmmdl_amber_logo correctly appends the full amber logo to the amber_script list."""
    session, uploadedFiles = mock_data
    amber_script_gen = AmberScriptGenerator(session, uploadedFiles)
    
    # Prepare an empty amber_script list
    amber_script = []

    # Call the method
    amber_script_gen.add_openmmdl_amber_logo(amber_script)

    # Define the expected Amber logo output
    expected_logo = """
#				      _              _               
#				     / \\   _ __ ___ | |__   ___ _ __ 
#				    / _ \\ | '_ ` _ \\| '_ \\ / _ \\ '__|
#				   / ___ \\| | | | | | |_) |  __/ |   
#				  /_/   \\_\\_| |_| |_|_.__/ \\___|_|    
        """

    # Verify that the full expected logo is appended
    assert len(amber_script) == 1
    assert amber_script[0] == expected_logo
    
    
def test_add_prot_receptor_type(base_mock_data):
    """Test if add_receptor_type correctly appends commands for protein receptor."""
    data = base_mock_data["protRcp"]
    amber_script_gen = AmberScriptGenerator(data["session"], data["uploadedFiles"])
    
    amber_script = []
    amber_script_gen.add_receptor_type(amber_script)

    expected_output = [
        "#!/bin/bash\n",
        "################################## Receptor ######################################\n",
        "rcp_nm=protein # the file name of ligand without suffix `pdb`",
        "rcp_ff=ff14SB",
        "\n"
    ]
    assert amber_script == expected_output

def test_add_dna_receptor_type(base_mock_data):
    """Test if add_receptor_type correctly appends commands for DNA receptor."""
    data = base_mock_data["dnaRcp"]
    amber_script_gen = AmberScriptGenerator(data["session"], data["uploadedFiles"])
    
    amber_script = []
    amber_script_gen.add_receptor_type(amber_script)

    expected_output = [
        "#!/bin/bash\n",
        "################################## Receptor ######################################\n",
        "rcp_nm=dna # the file name of ligand without suffix `pdb`",
        "rcp_ff=bsc1",
        "\n"
    ]
    assert amber_script == expected_output

def test_add_rna_receptor_type(base_mock_data):
    """Test if add_receptor_type correctly appends commands for RNA receptor."""
    data = base_mock_data["rnaRcp"]
    amber_script_gen = AmberScriptGenerator(data["session"], data["uploadedFiles"])
    
    amber_script = []
    amber_script_gen.add_receptor_type(amber_script)

    expected_output = [
        "#!/bin/bash\n",
        "################################## Receptor ######################################\n",
        "rcp_nm=rna # the file name of ligand without suffix `pdb`",
        "rcp_ff=ff99SB",
        "\n"
    ]
    assert amber_script == expected_output

def test_add_carbo_receptor_type(base_mock_data):
    """Test if add_receptor_type correctly appends commands for carbohydrate receptor."""
    data = base_mock_data["carboRcp"]
    amber_script_gen = AmberScriptGenerator(data["session"], data["uploadedFiles"])
    
    amber_script = []
    amber_script_gen.add_receptor_type(amber_script)

    expected_output = [
        "#!/bin/bash\n",
        "################################## Receptor ######################################\n",
        "rcp_nm=carbo # the file name of ligand without suffix `pdb`",
        "rcp_ff=GLYCAM",
        "\n"
    ]
    assert amber_script == expected_output
    
def test_add_clean_pdb_commands(mock_data):
    """Test if add_clean_pdb_commands correctly appends commands to clean the PDB file."""
    session, uploadedFiles = mock_data["session"], mock_data["uploadedFiles"]
    amber_script_gen = AmberScriptGenerator(session, uploadedFiles)

    amber_script = []
    amber_script_gen.add_clean_pdb_commands(amber_script)

    expected_output = [
        "## Clean the PDB file by pdb4amber",
        "pdb4amber -i ${rcp_nm}.pdb -o ${rcp_nm}_amber.pdb",
        """
## `tleap` requires that all residues and atoms have appropriate types to ensure compatibility with the specified force field.
## To avoid `tleap` failing, we delete non-essential atoms, such as hydrogens, but preserve important atoms like carbon and nitrogen within the caps residues.
## Don' worry about the missing atoms as tleap has the capability to reconstruct them automatically.""",
        """awk '! ($2 ~ "(CH3|HH31|HH32|HH33)" || $3 ~ "(CH3|HH31|HH32|HH33)" )' ${rcp_nm}_amber.pdb > ${rcp_nm}_amber_f.pdb""",
        "grep -v '^CONECT' ${rcp_nm}_amber_f.pdb > ${rcp_nm}_cnt_rmv.pdb\n"
    ]

    assert amber_script == expected_output


def test_add_ligand_commands_nmLig(mock_data):
    """Test if add_ligand_commands correctly appends commands for a normal ligand (nmLig)."""
    session = {
        "nmLig": True,
        "charge_method": "bcc",
        "charge_value": "-1",
        "lig_ff": "gaff"
    }
    uploadedFiles = {
        "nmLigFile": [("file5", "ligand.pdb")]
    }
    amber_script_gen = AmberScriptGenerator(session, uploadedFiles)

    amber_script = []
    amber_script_gen.add_ligand_commands(amber_script)

    expected_output = [
        "################################## Ligand ######################################",
        "# Normal Ligand that is compatible with GAFF force field",
        "nmLigFile=ligand # the file name of ligand without suffix `.pdb` or `.sdf`",
        "obabel ${nmLigFile}.pdb -O ${nmLigFile}.sdf -p # convert to sdf file for openmmdl_analysis, -p: add hydrogens appropriate for pH7.4",
        "charge_method=bcc # refers to the charge method that antechamber will adopt",
        "charge_value=-1 # Enter the net molecular charge of the ligand as integer (e.g. 1 or -2)",
        "lig_ff=gaff # Ligand force field\n",
        "## Clean the PDB file by pdb4amber",
        "pdb4amber -i ${nmLigFile}.pdb -o ${nmLigFile}_amber.pdb\n",
        "## Generate a prepc file and an additional frcmod file by `antechamber`",
        "antechamber -fi pdb -fo prepc -i ${nmLigFile}_amber.pdb -o ${nmLigFile}.prepc -c ${charge_method} -at ${lig_ff} -nc ${charge_value} -pf y",
        "parmchk2 -f prepc -i ${nmLigFile}.prepc -o ${nmLigFile}.frcmod\n",
        "## Rename ligand pdb",
        "antechamber -i ${nmLigFile}.prepc -fi prepc -o rename_${nmLigFile}.pdb -fo pdb\n"
    ]

    assert amber_script == expected_output


def test_add_ligand_commands_spLig(mock_data):
    """Test if add_ligand_commands correctly appends commands for a special ligand (spLig)."""
    session = {
        "spLig": True
    }
    uploadedFiles = {
        "spLigFile": [("file6", "ligand.pdb")],
        "prepcFile": [("file7", "ligand.prepc")],
        "frcmodFile": [("file8", "ligand.frcmod")]
    }
    amber_script_gen = AmberScriptGenerator(session, uploadedFiles)

    amber_script = []
    amber_script_gen.add_ligand_commands(amber_script)

    expected_output = [
        "################################## Ligand ######################################",
        "# Special Ligand that is incompatible with GAFF force field",
        "spLigFile=ligand # the file name of ligand without suffix `.pdb`",
        "prepc=ligand # the file name without suffix `prepc`",
        "frcmod=ligand # the file name without suffix `frcmod`\n",
        "## Clean the PDB file by pdb4amber",
        "pdb4amber -i ${spLigFile}.pdb -o ${spLigFile}_amber.pdb\n",
        "spLigName=$(awk 'NR==1 {print $4}' ${spLigFile}_amber.pdb)\n"
    ]

    assert amber_script == expected_output

def test_add_combine_components_commands(mock_data):
    """Test if add_combine_components_commands correctly appends commands to combine all components."""
    session, uploadedFiles = mock_data["session"], mock_data["uploadedFiles"]
    amber_script_gen = AmberScriptGenerator(session, uploadedFiles)

    amber_script = []
    amber_script_gen.add_combine_components_commands(amber_script)

    expected_output = [
        "######################  Combine All Components to Be Modelled ####################",
        "cat > tleap.combine.in <<EOF\n",
        "source ${rcp_ff}",
        "source leaprc.${lig_ff}",
        "loadamberprep ${nmLigFile}.prepc",
        "loadamberparams ${nmLigFile}.frcmod\n",
        "loadamberprep ${prepc}.prepc",
        "loadamberparams ${frcmod}.frcmod\n",
        "rcp = loadpdb ${rcp_nm}_cnt_rmv.pdb",
        "nmLig = loadpdb rename_${nmLigFile}.pdb ",
        "spLig = loadpdb ${spLigFile}_amber.pdb ",
        "comp = combine{rcp nmLig spLig}",
        "savepdb comp comp.pdb",
        "quit\nEOF\n",
        "tleap -s -f tleap.combine.in > tleap.combine.out",
        "grep -v '^CONECT' comp.pdb > comp_cnt_rmv.pdb\n"
    ]

    assert amber_script == expected_output


def test_add_solvation_commands(mock_data):
    """Test if add_solvation_commands correctly appends commands for solvation settings."""
    session = {
        "addType": "addWater",
        "boxType": "cube",
        "dist": "10.0"
    }
    uploadedFiles = {}
    amber_script_gen = AmberScriptGenerator(session, uploadedFiles)

    amber_script = []
    amber_script_gen.add_solvation_commands(amber_script)

    expected_output = [
        "boxType=solvatebox # `solvatebox`, a command in tleap, creates a cubic box ",
        "dist=10.0 # the minimum distance between any atom originally present in solute and the edge of the periodic box."
    ]

    assert amber_script == expected_output
