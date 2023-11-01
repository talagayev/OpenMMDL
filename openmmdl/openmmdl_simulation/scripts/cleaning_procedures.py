import os
import shutil

def cleanup(protein_name):
    """
    Cleans up the PDB Reporter Output File and MDTraj Files of the performed simulation

    Parameters
    ----------
    protein_name: str
        Name of the protein pdb.

    Returns
    ----------
    None
    """ 
    print("Cleaning Up :)")
    try:
        os.remove(f'output_{protein_name}')
        os.remove(f'centered_old_coordinates.pdb')
        os.remove(f'centered_old_coordinates.dcd')
    except FileNotFoundError:
        print("One or more files not found. Cleanup skipped.")
    print("Cleanup is done.")
    
def create_directory_if_not_exists(directory_path):
    """
    Create a directory if it doesn't exist, or overwrite it if already does.

    Parameters
    ----------
    directory_path: str
        Path of the directory that you want to create.

    Returns
    ----------
    None
    """ 
    if not os.path.exists(directory_path):
        os.mkdir(directory_path)
    else:
        shutil.rmtree(directory_path)
        os.mkdir(directory_path)

def copy_file(src, dest):
    """
    Copy a file to the destination path.

    Parameters
    ----------
    src: str
        Path of the file that needs to be copied.
    dest: str
        Path of destination where the file needs to be copied to.

    Returns
    ----------
    None
    """ 
    if os.path.exists(src):
        shutil.copy(src, dest)

def organize_files(source, destination):
    """
    Organizes the files and moves them from the source to the destination directory.

    Parameters
    ----------
    source: str
        Path of the file that needs to be moved.
    destination: str
        Path of destination where the file needs to be moved to.

    Returns
    ----------
    None
    """ 
    for file in source:
        if os.path.exists(file):
            os.rename(file, os.path.join(destination, os.path.basename(file)))

def post_md_file_movement(protein_name, prmtop=None, inpcrd=None, ligand=None):
    """
    Organizes and moves the files after the MD simulation to their respective directories.
    Parameters
    ----------
    protein_name : str
        Name of the protein pdb.
    prmtop : str (optional)
        Path to the AMBER topology file.
    inpcrd : str (optional)
        Path to the AMBER coordinate file.
    ligand : str (optional)
        Path to the ligand file.

    Returns
    ----------
    None
    """ 
    # Create necessary directories
# Create necessary directories
    create_directory_if_not_exists("Input_Files")
    create_directory_if_not_exists("MD_Files")
    create_directory_if_not_exists("MD_Files/Pre_MD")
    create_directory_if_not_exists("MD_Files/Minimization_Equilibration")
    create_directory_if_not_exists("MD_Files/MD_Output")
    create_directory_if_not_exists("MD_Postprocessing")
    create_directory_if_not_exists("Final_Output")
    create_directory_if_not_exists("Final_Output/All_Atoms")
    create_directory_if_not_exists("Final_Output/Prot_Lig")
    create_directory_if_not_exists("Checkpoints")

    # Print if directories have been created
    for directory in ["Input_Files", "MD_Files", "MD_Files/Pre_MD", "MD_Files/Minimization_Equilibration",
                     "MD_Files/MD_Output", "MD_Postprocessing", "Final_Output", "Final_Output/All_Atoms",
                     "Final_Output/Prot_Lig", "Checkpoints"]:
        if os.path.exists(directory):
            print(f"Directory '{directory}' exists.")
        else:
            print(f"Directory '{directory}' does not exist.")

    # Move input files
    if ligand:
        copy_file(ligand, "Final_Output/All_Atoms")
        copy_file(ligand, "Final_Output/Prot_Lig")
        copy_file(ligand, "Input_Files")  # Add this line to copy the ligand to Input_Files
        print("Copied ligand to Final_Output/All_Atoms, Final_Output/Prot_Lig, and Input_Files")

        ligand_dir = "Final_Output/All_Atoms"
        if os.path.exists(os.path.join(ligand_dir, os.path.basename(ligand))):
            print(f"File '{os.path.basename(ligand)}' exists in '{ligand_dir}'")
        else:
            print(f"File '{os.path.basename(ligand)}' does not exist in '{ligand_dir}'")

        input_files_dir = "Input_Files"
        if os.path.exists(os.path.join(input_files_dir, os.path.basename(ligand))):
            print(f"File '{os.path.basename(ligand)}' exists in '{input_files_dir}'")
        else:
            print(f"File '{os.path.basename(ligand)}' does not exist in '{input_files_dir}'")

    if protein_name:
        copy_file(protein_name, "Input_Files")
        print(f"Copied {os.path.basename(protein_name)} to Input_Files")

        input_files_dir = "Input_Files"
        if os.path.exists(os.path.join(input_files_dir, os.path.basename(protein_name))):
            print(f"File '{os.path.basename(protein_name)}' exists in '{input_files_dir}'")
        else:
            print(f"File '{os.path.basename(protein_name)}' does not exist in '{input_files_dir}'")

    if prmtop:
        copy_file(prmtop, "Input_Files")
        print(f"Copied {os.path.basename(prmtop)} to Input_Files")

        input_files_dir = "Input_Files"
        if os.path.exists(os.path.join(input_files_dir, os.path.basename(prmtop))):
            print(f"File '{os.path.basename(prmtop)}' exists in '{input_files_dir}'")
        else:
            print(f"File '{os.path.basename(prmtop)}' does not exist in '{input_files_dir}'")

    if inpcrd:
        copy_file(inpcrd, "Input_Files")
        print(f"Copied {os.path.basename(inpcrd)} to Input_Files")

        input_files_dir = "Input_Files"
        if os.path.exists(os.path.join(input_files_dir, os.path.basename(inpcrd))):
            print(f"File '{os.path.basename(inpcrd)}' exists in '{input_files_dir}'")
        else:
            print(f"File '{os.path.basename(inpcrd)}' does not exist in '{input_files_dir}'")

    # Organize pre-MD files
    source_pre_md_files = ["prepared_no_solvent_", "solvent_padding_", "solvent_absolute_", "membrane_"]
    destination_pre_md = "MD_Files/Pre_MD"
    organize_files([f"{prefix}{protein_name}" for prefix in source_pre_md_files], destination_pre_md)

    # Organize topology files after minimization and equilibration
    source_topology_files = ["Energyminimization_", "Equilibration_"]
    destination_topology = "MD_Files/Minimization_Equilibration"
    organize_files([f"{prefix}{protein_name}" for prefix in source_topology_files], destination_topology)

    # Organize simulation output files
    organize_files([f"output_{protein_name}", "trajectory.dcd"], "MD_Files/MD_Output")

    # Organize MDtraj and MDAnalysis files
    organize_files(["centered_old_coordinates_top.pdb", "centered_old_coordinates.dcd", "centered_old_coordinates_top.gro", "centered_old_coordinates.xtc"], "MD_Postprocessing")
    organize_files(["centered_top.pdb", "centered_traj.dcd", "centered_top.gro", "centered_traj.xtc"], "Final_Output/All_Atoms")
    organize_files(["prot_lig_top.pdb", "prot_lig_traj.dcd", "prot_lig_top.gro", "prot_lig_traj.xtc"], "Final_Output/Prot_Lig")

    # Organize checkpoint files
    organize_files(["checkpoint.chk", "10x_checkpoint.chk", "100x_checkpoint.chk"], "Checkpoints")
