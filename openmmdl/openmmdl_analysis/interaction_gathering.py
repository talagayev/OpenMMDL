import os
import pandas as pd
import MDAnalysis as mda
from tqdm import tqdm
from plip.basic import config
from plip.structure.preparation import PDBComplex, LigandFinder, Mol, PLInteraction
from plip.exchange.report import BindingSiteReport
from multiprocessing import Pool
from functools import partial

config.KEEPMOD = True


def characterize_complex(pdb_file: str, binding_site_id: str) -> PLInteraction:
    """Characterize the protein-ligand complex and return their interaction set

    Args:
        pdb_file (str): A string, which represents the path to the PDB File
        binding_site_id (str): A string that specifies the identifier of the binding site

    Returns:
        PLInteraction: A object representing the interactions if. If Binding site is not found returns None
    """
    pdb_complex = PDBComplex()
    pdb_complex.load_pdb(pdb_file)
    for ligand in pdb_complex.ligands:
        if (
            ":".join([ligand.hetid, ligand.chain, str(ligand.position)])
            == binding_site_id
        ):
            pdb_complex.characterize_complex(ligand)

    return pdb_complex.interaction_sets[binding_site_id]


def retrieve_plip_interactions(pdb_file, lig_name):
    """Retrieves the interactions from PLIP.

    Args:
        pdb_file (str): The path of the PDB file of the complex.
        lig_name (str): Name of the Ligand in the complex topology that will be analyzed.

    Returns:
        dict: A dictionary of the binding sites and the interactions.
    """
    protlig = PDBComplex()
    protlig.load_pdb(pdb_file)  # load the pdb file
    for ligand in protlig.ligands:
        if str(ligand.longname) == lig_name:
            protlig.characterize_complex(
                ligand
            )  # find ligands and analyze interactions
    sites = {}
    # loop over binding sites
    for key, site in sorted(protlig.interaction_sets.items()):
        binding_site = BindingSiteReport(site)  # collect data about interactions
        # tuples of *_features and *_info will be converted to pandas DataFrame
        keys = (
            "hydrophobic",
            "hbond",
            "waterbridge",
            "saltbridge",
            "pistacking",
            "pication",
            "halogen",
            "metal",
        )
        # interactions is a dictionary which contains relevant information for each
        # of the possible interactions: hydrophobic, hbond, etc. in the considered
        # binding site.
        interactions = {
            k: [getattr(binding_site, k + "_features")]
            + getattr(binding_site, k + "_info")
            for k in keys
        }
        sites[key] = interactions

    return sites


def retrieve_plip_interactions_peptide(pdb_file, peptide):
    """Retrives the interactions from PLIP for a peptide.

    Args:
        pdb_file (str): The path of the PDB file of the complex.
        peptide (str): Chainid of the peptide that will be analyzed.

    Returns:
        dict: A dictionary of the binding sites and the interactions.
    """
    protlig = PDBComplex()
    protlig.load_pdb(pdb_file)  # load the pdb file
    protlig.characterize_complex(
        protlig.ligands[-1]
    )  # find ligands and analyze interactions
    sites = {}
    # loop over binding sites
    for key, site in sorted(protlig.interaction_sets.items()):
        binding_site = BindingSiteReport(site)  # collect data about interactions
        # tuples of *_features and *_info will be converted to pandas DataFrame
        keys = (
            "hydrophobic",
            "hbond",
            "waterbridge",
            "saltbridge",
            "pistacking",
            "pication",
            "halogen",
            "metal",
        )
        # interactions is a dictionary which contains relevant information for each
        # of the possible interactions: hydrophobic, hbond, etc. in the considered
        # binding site.
        interactions = {
            k: [getattr(binding_site, k + "_features")]
            + getattr(binding_site, k + "_info")
            for k in keys
        }
        sites[key] = interactions

    return sites


def create_df_from_binding_site(selected_site_interactions, interaction_type="hbond"):
    """Creates a data frame from a binding site and interaction type.

    Args:
        selected_site_interactions (dict): Precaluclated interactions from PLIP for the selected site
        interaction_type (str, optional): The interaction type of interest (default set to hydrogen bond). Defaults to "hbond".

    Returns:
        pandas dataframe: DataFrame with information retrieved from PLIP.
    """
    # check if interaction type is valid:
    valid_types = [
        "hydrophobic",
        "hbond",
        "waterbridge",
        "saltbridge",
        "pistacking",
        "pication",
        "halogen",
        "metal",
    ]

    if interaction_type not in valid_types:
        print(
            "\033[1m!!! Wrong interaction type specified. Hbond is chosen by default!!!\033[0m\n"
        )
        interaction_type = "hbond"

    df = pd.DataFrame.from_records(
        # data is stored AFTER the column names
        selected_site_interactions[interaction_type][1:],
        # column names are always the first element
        columns=selected_site_interactions[interaction_type][0],
    )
    return df


def change_lig_to_residue(file_path, old_residue_name, new_residue_name):
    """Reformats the topology file to change the ligand to a residue. This is needed for interactions with special ligands such as metal ions.

    Args:
        file_path (str): Filepath of the topology file.
        old_residue_name (str): Residue name of the ligand.
        new_residue_name (str): New residue name of the ligand now changed to mimic a residue.
    """
    with open(file_path, "r") as file:
        lines = file.readlines()

    with open(file_path, "w") as file:
        for line in lines:
            if line.startswith("HETATM") or line.startswith("ATOM"):
                # Assuming the standard PDB format for simplicity
                # You may need to adapt this part based on your specific PDB file
                atom_name = line[12:16].strip()
                residue_name = line[17:20].strip()

                # Check if the residue name matches the one to be changed
                if residue_name == old_residue_name:
                    # Change the residue name to the new one
                    modified_line = line[:17] + new_residue_name + line[20:]
                    file.write(modified_line)
                else:
                    file.write(line)
            else:
                file.write(line)


def process_frame(frame, pdb_md, lig_name, special=None, peptide=None):
    """Process a single frame of MD simulation.

    Args:
        frame (int): The number of the frame that will be processed.
        pdb_md (mda universe): The MDAnalysis Universe class representation of the topology and the trajectory of the file that is being processed.
        lig_name (str): Name of the ligand in the complex that will be analyzed.
        special (str, optional): Name of the special ligand in the complex that will be analyzed. Defaults to None.
        peptide (srt, optional): Chainid of the peptide that will be analyzed. Defaults to None.

    Returns:
        pandas dataframe: A dataframe conatining the interaction data for the processed frame.
    """
    atoms_selected = pdb_md.select_atoms(
        f"protein or nucleic or resname {lig_name} or (resname HOH and around 10 resname {lig_name}) or resname {special}"
    )
    for num in pdb_md.trajectory[(frame) : (frame + 1)]:
        atoms_selected.write(f"processing_frame_{frame}.pdb")
    if peptide is None:
        interactions_by_site = retrieve_plip_interactions(
            f"processing_frame_{frame}.pdb", lig_name
        )
        index_of_selected_site = -1
        selected_site = list(interactions_by_site.keys())[index_of_selected_site]

        interaction_types = [
            "hydrophobic",
            "hbond",
            "waterbridge",
            "saltbridge",
            "pistacking",
            "pication",
            "halogen",
            "metal",
        ]

        interaction_list = pd.DataFrame()
        for interaction_type in interaction_types:
            tmp_interaction = create_df_from_binding_site(
                interactions_by_site[selected_site], interaction_type=interaction_type
            )
            tmp_interaction["FRAME"] = int(frame)
            tmp_interaction["INTERACTION"] = interaction_type
            interaction_list = pd.concat([interaction_list, tmp_interaction])
        if os.path.exists(f"processing_frame_{frame}.pdb"):
            os.remove(f"processing_frame_{frame}.pdb")

        if special is not None:
            combi_lig_special = mda.Universe("ligand_special.pdb")
            complex = mda.Universe("complex.pdb")
            complex_all = complex.select_atoms("all")
            result = process_frame_special(frame, pdb_md, lig_name, special)
            results_df = pd.concat(result, ignore_index=True)
            results_df = results_df[results_df["LOCATION"] == "protein.sidechain"]
            results_df["RESTYPE"] = results_df["RESTYPE"].replace(
                ["HIS", "SER", "CYS"], lig_name
            )
            results_df["LOCATION"] = results_df["LOCATION"].replace(
                "protein.sidechain", "ligand"
            )
            updated_target_idx = []

            for index, row in results_df.iterrows():
                ligand_special_int_nr = int(row["TARGET_IDX"])
                ligand_special_int_nr_atom = combi_lig_special.select_atoms(
                    f"id {ligand_special_int_nr}"
                )
                for atom in ligand_special_int_nr_atom:
                    atom_name = atom.name
                    # Adjust atom_name based on the specified conditions
                    if atom_name in ["N", "C", "O", "S"]:
                        atom_name = f"{atom_name}1"
                    else:
                        # Assuming the format is a single letter followed by a number
                        base_name, atom_number = atom_name[:-1], int(atom_name[-1])
                        new_atom_number = atom_number + 1
                        atom_name = f"{base_name}{new_atom_number}"
                    for complex_atom in complex_all:
                        complex_atom_name = complex_atom.name
                        if atom_name == complex_atom_name:
                            true_number = complex_atom.id
                            break  # Exit the loop once a match is found
                    updated_target_idx.append(true_number)

            # Update 'TARGET_IDX' in interaction_list
            results_df["TARGET_IDX"] = updated_target_idx
            interaction_list["TARGET_IDX"] = interaction_list["TARGET_IDX"]

            # Concatenate the updated results_df to interaction_list
            interaction_list = pd.concat([interaction_list, results_df])
    if peptide is not None:
        interactions_by_site = retrieve_plip_interactions_peptide(
            f"processing_frame_{frame}.pdb", peptide
        )
        index_of_selected_site = -1
        selected_site = list(interactions_by_site.keys())[index_of_selected_site]

        interaction_types = [
            "hydrophobic",
            "hbond",
            "waterbridge",
            "saltbridge",
            "pistacking",
            "pication",
            "halogen",
            "metal",
        ]

        interaction_list = pd.DataFrame()
        for interaction_type in interaction_types:
            tmp_interaction = create_df_from_binding_site(
                interactions_by_site[selected_site], interaction_type=interaction_type
            )
            tmp_interaction["FRAME"] = int(frame)
            tmp_interaction["INTERACTION"] = interaction_type
            interaction_list = pd.concat([interaction_list, tmp_interaction])
        if os.path.exists(f"processing_frame_{frame}.pdb"):
            os.remove(f"processing_frame_{frame}.pdb")

    return interaction_list


def process_frame_special(frame, pdb_md, lig_name, special=None):
    """Function extension of process_frame to process special ligands.

    Args:
        frame (int): Number of the frame that will be processed.
        pdb_md (mda universe): MDA Universe class representation of the topology and the trajectory of the file that is being processed.
        lig_name (str): Name of the ligand in the complex that will be analyzed.
        special (str, optional): Name of the special ligand that will be analysed. Defaults to None.

    Returns:
        list: list of dataframes containing the interaction data for the processed frame with the special ligand.
    """
    res_renaming = ["HIS", "SER", "CYS"]
    interaction_dfs = []
    for res in res_renaming:
        pdb_md.trajectory[frame]
        atoms_selected = pdb_md.select_atoms(f"resname {lig_name} or resname {special}")
        atoms_selected.write(f"processing_frame_{frame}.pdb")
        change_lig_to_residue(f"processing_frame_{frame}.pdb", lig_name, res)
        interactions_by_site = retrieve_plip_interactions(
            f"processing_frame_{frame}.pdb", special
        )
        index_of_selected_site = -1
        selected_site = list(interactions_by_site.keys())[index_of_selected_site]
        interaction_types = ["metal"]
        interaction_list = pd.DataFrame()
        for interaction_type in interaction_types:
            tmp_interaction = create_df_from_binding_site(
                interactions_by_site[selected_site], interaction_type=interaction_type
            )
            tmp_interaction["FRAME"] = int(frame)
            tmp_interaction["INTERACTION"] = interaction_type
            interaction_list = pd.concat([interaction_list, tmp_interaction])
        interaction_dfs.append(interaction_list)
        os.remove(f"processing_frame_{frame}.pdb")
    return interaction_dfs


def process_frame_wrapper(args):
    """Wrapper for the MD Trajectory procession.

    Args:
        args (tuple): Tuple containing (frame_idx: int - number of the frame to be processed, pdb_md: mda universe - MDA Universe class representation of the topology and the trajectory of the file that is being processed, lig_name: str - Name of the ligand in the complex that will be analyzed, special_ligand: str - Name of the special ligand that will be analysed, peptide: str - Chainid of the peptide that will be analyzed)

    Returns:
        tuple: tuple containing the frame index and the result of from the process_frame function.
    """
    frame_idx, pdb_md, lig_name, special_ligand, peptide = args

    return frame_idx, process_frame(
        frame_idx, pdb_md, lig_name, special_ligand, peptide
    )


def process_trajectory(
    pdb_md, dataframe, num_processes, lig_name, special_ligand, peptide
):
    """Process protein-ligand trajectory with multiple CPUs in parallel.

    Args:
        pdb_md (mda universe): MDAnalysis Universe class representation of the topology and the trajectory of the file that is being processed.
        dataframe (pandas dataframe): Name of a CSV file as str, where the interaction data will be read from if not None.
        num_processes (int): The number of CPUs that will be used for the processing of the protein-ligand trajectory. Defaults to half of the CPUs in the system.
        lig_name (str): Name of the Ligand in the complex that will be analyzed.
        special_ligand (str): Name of the special ligand in the complex that will be analyzed.
        peptide (str): Chainid of the peptide that will be analyzed.

    Returns:
        pandas dataframe: A DataFrame containing all the protein-ligand interaction data from the whole trajectory.
    """
    if dataframe is None:
        print("\033[1mProcessing protein-ligand trajectory\033[0m")
        print(f"\033[1mUsing {num_processes} CPUs\033[0m")
        total_frames = len(pdb_md.trajectory) - 1

        with Pool(processes=num_processes) as pool:
            frame_args = [
                (i, pdb_md, lig_name, special_ligand, peptide)
                for i in range(1, total_frames + 1)
            ]

            # Initialize the progress bar with the total number of frames
            pbar = tqdm(total=total_frames, ascii=True, desc="Analyzing frames")

            results = []
            for result in pool.imap(process_frame_wrapper, frame_args):
                results.append(result)
                pbar.update(1)  # Update the progress manually

        # Close the progress bar
        pbar.close()

        # Extract the results and sort them by frame index
        results.sort(key=lambda x: x[0])
        interaction_lists = [result[1] for result in results]

        interaction_list = pd.concat(interaction_lists)

        interaction_list.to_csv("interactions_gathered.csv")

    elif dataframe is not None:
        print(f"\033[1mGathering data from {dataframe}\033[0m")
        interaction_tmp = pd.read_csv(dataframe)
        interaction_list = interaction_tmp.drop(interaction_tmp.columns[0], axis=1)

    print("\033[1mProtein-ligand trajectory processed\033[0m")

    return interaction_list


def fill_missing_frames(df, md_len):
    """Fills the frames with no interactions in the DataFrame with placeholder values.

    Args:
        df (pandas dataframe): The input DataFrame with frames that have no Interactions
        md_len (int): The value that indicates the number of frames, thus allowing the function to loop through the DataFrame

    Returns:
        pandas dataframe: DataFrame with placeholder values in the frames with no interactions.
    """
    # Create a set containing all unique values in the 'FRAME' column
    existing_frames = set(df["FRAME"])

    # Create a list to store new rows for missing numbers
    missing_rows = []

    # Iterate through numbers from 0 to md_len
    for frame_number in range(1, md_len):
        if frame_number not in existing_frames:
            # Create a new row with 'FRAME' set to the missing number and other columns set to "skip"
            missing_row = {"FRAME": frame_number}
            for col in df.columns:
                if col != "FRAME":
                    missing_row[col] = "skip"
            missing_rows.append(missing_row)

    # Concatenate the missing rows with the original DataFrame
    df = pd.concat([df, pd.DataFrame(missing_rows)], ignore_index=True)

    # Sort the DataFrame by the 'FRAME' column
    df.sort_values(by="FRAME", inplace=True)

    return df
