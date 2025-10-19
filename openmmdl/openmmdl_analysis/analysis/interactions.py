import re
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple

import MDAnalysis as mda
import pandas as pd
from multiprocessing import Pool
from prolif import Fingerprint
from tqdm import tqdm


@dataclass
class _AtomInfo:
    """Container holding atom level information extracted from ProLIF."""

    index: Optional[int]
    resid: Optional[int]
    resname: Optional[str]
    chain: Optional[str]
    name: Optional[str]


def _get_atom_info(atom) -> _AtomInfo:
    """Return a normalized :class:`_AtomInfo` instance for a ProLIF atom entry."""

    if atom is None:
        return _AtomInfo(None, None, None, None, None)

    if hasattr(atom, "resid"):
        resid = getattr(atom, "resid", None)
        resname = getattr(atom, "resname", None)
        chain = getattr(atom, "chainID", getattr(atom, "segid", None))
        name = getattr(atom, "name", None)
        index = getattr(atom, "index", getattr(atom, "id", None))
        if index is not None:
            index = int(index)
        if resid is not None:
            resid = int(resid)
        return _AtomInfo(index, resid, resname, chain, name)

    if isinstance(atom, dict):
        index = atom.get("index") or atom.get("id") or atom.get("serial")
        resid = atom.get("resid") or atom.get("resnum")
        resname = atom.get("resname") or atom.get("restype")
        chain = atom.get("chain") or atom.get("chain_id")
        name = atom.get("name") or atom.get("atomname")
        if index is not None:
            index = int(index)
        if resid is not None:
            resid = int(resid)
        return _AtomInfo(index, resid, resname, chain, name)

    if isinstance(atom, (tuple, list)) and atom:
        return _get_atom_info(atom[0])

    return _AtomInfo(None, None, None, None, None)


def _parse_residue_label(label: str) -> Tuple[Optional[str], Optional[int], Optional[str]]:
    """Parse a residue label from ProLIF into chain, resid and resname."""

    if not label:
        return None, None, None

    # Expected formats include "A:LEU132" or "LEU132" or "A:LEU:132" or "LIG1"
    # Use regex to capture the residue name and number.
    match = re.search(
        r"(?:(?P<chain>[A-Za-z0-9]))?:?(?P<resname>[A-Za-z0-9]+)(?:[^0-9]*(?P<resid>-?\d+))?",
        str(label),
    )
    if match:
        chain = match.group("chain")
        resid = match.group("resid")
        resname = match.group("resname")
        return chain, int(resid) if resid else None, resname

    return None, None, str(label)


def _normalise_match_entry(entry) -> Tuple[List[_AtomInfo], List[_AtomInfo], Dict[str, object]]:
    """Return normalized lists of protein/ligand atom data and metadata for a match."""

    if isinstance(entry, dict):
        protein_atoms = entry.get("protein_atoms") or entry.get("protein") or []
        ligand_atoms = entry.get("ligand_atoms") or entry.get("ligand") or []
        metadata = entry.get("metadata") or entry.get("info") or {}
    elif isinstance(entry, (tuple, list)):
        protein_atoms = entry[0] if len(entry) > 0 else []
        ligand_atoms = entry[1] if len(entry) > 1 else []
        metadata = {}
    else:
        protein_atoms = []
        ligand_atoms = []
        metadata = {}

    protein_atoms = [
        _get_atom_info(atom)
        for atom in (protein_atoms if isinstance(protein_atoms, (list, tuple)) else [protein_atoms])
    ]
    ligand_atoms = [
        _get_atom_info(atom)
        for atom in (ligand_atoms if isinstance(ligand_atoms, (list, tuple)) else [ligand_atoms])
    ]

    return protein_atoms, ligand_atoms, metadata


def _create_empty_row() -> Dict[str, object]:
    """Return an empty interaction row initialised with expected columns."""

    return {
        "FRAME": None,
        "INTERACTION": None,
        "RESNR": None,
        "RESTYPE": None,
        "RESCHAIN": None,
        "RESNR_LIG": None,
        "RESTYPE_LIG": None,
        "RESCHAIN_LIG": None,
        "LOCATION": None,
        "Prot_partner": None,
        "LIGCARBONIDX": None,
        "PROTISDON": None,
        "ACCEPTORIDX": None,
        "DONORIDX": None,
        "ACCEPTOR_IDX": None,
        "DONOR_IDX": None,
        "LIG_IDX_LIST": None,
        "LIG_GROUP": None,
        "PROTISPOS": None,
        "TARGET_IDX": None,
        "METAL_TYPE": None,
        "COORDINATION": None,
        "DON_IDX": None,
        "DONORTYPE": None,
    }


def _format_atom_indices(atoms: Iterable[_AtomInfo]) -> Optional[object]:
    """Convert atom indices into a representation compatible with legacy consumers."""

    indices = [atom.index for atom in atoms if atom.index is not None]
    if not indices:
        return None
    unique = sorted(set(indices))
    if len(unique) == 1:
        return unique[0]
    return ",".join(str(idx) for idx in unique)



class InteractionAnalyzer:
    """
    Analyzes molecular interactions between a protein and a ligand/peptide
    throughout an MD trajectory using ProLIF.

    Attributes
    ----------
    pdb_md : mda.Universe
        MDAnalysis Universe object representing the topology and trajectory.
    dataframe : str or None
        Path to an existing interaction CSV file. If None, the trajectory will be processed anew.
    num_processes : int
        Number of CPU cores to use for parallel frame analysis.
    lig_name : str
        Residue name of the ligand in the complex.
    special_ligand : str
        Residue name for special ligands like metal ions (optional).
    peptide : str
        Chain ID of the peptide ligand (optional).
    md_len : int
        Number of frames in the trajectory.
    interaction_list : pd.DataFrame
        DataFrame storing the extracted interactions across the trajectory.
    """
    def __init__(
        self,
        pdb_md,
        dataframe,
        num_processes,
        lig_name,
        special_ligand,
        peptide,
        md_len,
    ):
        self.pdb_md = pdb_md
        self.dataframe = dataframe
        self.num_processes = num_processes
        self.lig_name = lig_name
        self.special = special_ligand
        self.peptide = peptide
        self.md_len = md_len
        self._prolif_interactions = [
            "Hydrophobic",
            "HBAcceptor",
            "HBDonor",
            "WaterBridge",
            "SaltBridge",
            "PiStacking",
            "PiCation",
            "CationPi",
            "HalogenBond",
            "MetalComplex",
        ]
        self.interaction_list = self._process_trajectory()

    def _select_ligand_group(self, frame: int) -> Tuple[mda.AtomGroup, mda.AtomGroup]:
        """Return ligand and receptor atom groups for the current frame."""

        self.pdb_md.trajectory[frame]
        if self.peptide is not None:
            ligand_selection = f"chainID {self.peptide}"
        else:
            ligand_selection = f"resname {self.lig_name}"

        base_selection = "protein or nucleic"
        if self.special is not None:
            base_selection += f" or resname {self.special}"

        receptor_selection = (
            f"({base_selection}) or (resname HOH and around 10 {ligand_selection})"
        )

        ligand_atoms = self.pdb_md.select_atoms(ligand_selection)
        receptor_atoms = self.pdb_md.select_atoms(receptor_selection)

        return ligand_atoms, receptor_atoms

    def _run_prolif(self, frame: int) -> pd.DataFrame:
        """Run ProLIF on a single frame and return the raw DataFrame."""

        ligand_atoms, receptor_atoms = self._select_ligand_group(frame)

        if len(ligand_atoms) == 0 or len(receptor_atoms) == 0:
            return pd.DataFrame()

        fingerprint = Fingerprint(interactions=self._prolif_interactions)

        # Create single-frame temporary universe to avoid sharing state across processes
        tmp_universe = mda.Merge(receptor_atoms, ligand_atoms)
        tmp_universe.load_new(self.pdb_md.trajectory.ts.copy())

        receptor_count = len(receptor_atoms)
        receptor_copy = tmp_universe.atoms[:receptor_count]
        ligand_copy = tmp_universe.atoms[receptor_count:]

        fingerprint.run(
            tmp_universe.trajectory[:],
            ligand=ligand_copy,
            protein=receptor_copy,
            progress=False,
        )

        try:
            prolif_df = fingerprint.to_dataframe(return_atoms=True)
        except TypeError:
            # Older versions of ProLIF use return_atoms keyword under ``keep``
            prolif_df = fingerprint.to_dataframe(keep=True)

        return prolif_df

    def _convert_prolif_dataframe(self, prolif_df: pd.DataFrame, frame: int) -> pd.DataFrame:
        """Convert ProLIF boolean dataframe into the legacy PLIP-like format."""

        if prolif_df.empty:
            return pd.DataFrame(columns=_create_empty_row().keys())

        rows: List[Dict[str, object]] = []
        current_frame = frame
        if isinstance(prolif_df.index, pd.MultiIndex):
            # When multiple frames are stored, select the requested frame
            if current_frame in prolif_df.index.get_level_values(0):
                frame_df = prolif_df.xs(current_frame, level=0)
            else:
                frame_df = prolif_df.iloc[0:1]
        else:
            frame_df = prolif_df

        if isinstance(frame_df, pd.Series):
            frame_iterable = frame_df.items()
        else:
            frame_iterable = frame_df.iloc[0].items()

        for column, value in frame_iterable:
            if value in (False, None) or (hasattr(value, "__len__") and len(value) == 0):
                continue

            if isinstance(value, (list, tuple)):
                match_entries = value
            else:
                match_entries = [value]

            if isinstance(column, tuple) and len(column) >= 3:
                interaction_label = column[0]
                protein_label = column[1]
                ligand_label = column[2]
            elif isinstance(column, tuple) and len(column) == 2:
                interaction_label = column[0]
                protein_label = column[1]
                ligand_label = None
            else:
                interaction_label = column
                protein_label = None
                ligand_label = None

            for entry in match_entries:
                protein_atoms, ligand_atoms, metadata = _normalise_match_entry(entry)

                row = _create_empty_row()
                row["FRAME"] = int(current_frame)

                chain, resid, resname = _parse_residue_label(protein_label)
                lig_chain, lig_resid, lig_resname = _parse_residue_label(ligand_label)

                if protein_atoms and protein_atoms[0].resid is not None:
                    resid = protein_atoms[0].resid
                if protein_atoms and protein_atoms[0].resname is not None:
                    resname = protein_atoms[0].resname
                if protein_atoms and protein_atoms[0].chain is not None:
                    chain = protein_atoms[0].chain

                if ligand_atoms and ligand_atoms[0].resid is not None:
                    lig_resid = ligand_atoms[0].resid
                if ligand_atoms and ligand_atoms[0].resname is not None:
                    lig_resname = ligand_atoms[0].resname
                if ligand_atoms and ligand_atoms[0].chain is not None:
                    lig_chain = ligand_atoms[0].chain

                row["RESNR"] = resid
                row["RESTYPE"] = resname
                row["RESCHAIN"] = chain
                row["RESNR_LIG"] = lig_resid
                row["RESTYPE_LIG"] = lig_resname
                row["RESCHAIN_LIG"] = lig_chain
                row["LOCATION"] = "protein.sidechain"

                if resid is not None and resname is not None and chain is not None:
                    row["Prot_partner"] = f"{resid}{resname}{chain}"

                interaction = str(interaction_label).lower()
                if interaction in ("hbdonor", "hbacceptor"):
                    interaction = "hbond"
                elif interaction == "metalcomplex":
                    interaction = "metal"
                elif interaction == "waterbridge":
                    interaction = "waterbridge"
                elif interaction == "halogenbond":
                    interaction = "halogen"
                elif interaction in ("pication", "cationpi"):
                    interaction = "pication"
                elif interaction == "pistacking":
                    interaction = "pistacking"
                elif interaction == "hydrophobic":
                    interaction = "hydrophobic"
                elif interaction == "saltbridge":
                    interaction = "saltbridge"

                row["INTERACTION"] = interaction

                if interaction == "hydrophobic":
                    row["LIGCARBONIDX"] = _format_atom_indices(ligand_atoms)
                elif interaction == "hbond":
                    prot_is_donor = str(interaction_label).lower() == "hbdonor"
                    row["PROTISDON"] = prot_is_donor
                    if prot_is_donor:
                        row["DONORIDX"] = _format_atom_indices(protein_atoms)
                        row["ACCEPTORIDX"] = _format_atom_indices(ligand_atoms)
                    else:
                        row["DONORIDX"] = _format_atom_indices(ligand_atoms)
                        row["ACCEPTORIDX"] = _format_atom_indices(protein_atoms)
                elif interaction == "waterbridge":
                    prot_is_donor = metadata.get("protein_is_donor")
                    if prot_is_donor is None:
                        prot_is_donor = str(interaction_label).lower() == "waterbridge"
                    row["PROTISDON"] = prot_is_donor
                    row["ACCEPTOR_IDX"] = _format_atom_indices(protein_atoms)
                    row["DONOR_IDX"] = _format_atom_indices(ligand_atoms)
                elif interaction == "pistacking":
                    row["LIG_IDX_LIST"] = _format_atom_indices(ligand_atoms)
                elif interaction == "pication":
                    row["LIG_IDX_LIST"] = _format_atom_indices(ligand_atoms)
                    lig_group = metadata.get("ligand_group")
                    if lig_group is None:
                        if str(interaction_label).lower() == "pication":
                            lig_group = "cation"
                        else:
                            lig_group = "pi"
                    row["LIG_GROUP"] = lig_group
                elif interaction == "saltbridge":
                    row["LIG_IDX_LIST"] = _format_atom_indices(ligand_atoms)
                    row["LIG_GROUP"] = metadata.get("ligand_group")
                    prot_charge = metadata.get("protein_charge")
                    if prot_charge is None and protein_atoms:
                        prot_charge = metadata.get("protein_type")
                    if isinstance(prot_charge, str):
                        prot_charge = prot_charge.lower() in ("positive", "pos", "+")
                    row["PROTISPOS"] = prot_charge
                elif interaction == "halogen":
                    row["DON_IDX"] = _format_atom_indices(ligand_atoms)
                    if ligand_atoms and ligand_atoms[0].name:
                        row["DONORTYPE"] = ligand_atoms[0].name
                elif interaction == "metal":
                    row["TARGET_IDX"] = _format_atom_indices(protein_atoms)
                    row["METAL_TYPE"] = metadata.get("metal_type") or (
                        ligand_atoms[0].resname if ligand_atoms else None
                    )
                    row["COORDINATION"] = metadata.get("coordination")
                    if self.special is not None:
                        row["RESTYPE_LIG"] = self.special

                rows.append(row)

        if rows:
            return pd.DataFrame(rows)

        return pd.DataFrame(columns=_create_empty_row().keys())

    def _process_frame(self, frame):
        """
        Process a single frame of MD simulation.

        Parameters
        ----------
        frame : int
            The number of the frame that will be processed.

        Returns
        -------
        pd.DataFrame
            A dataframe conatining the interaction data for the processed frame.
        """
        prolif_df = self._run_prolif(frame)
        interaction_list = self._convert_prolif_dataframe(prolif_df, frame)
        return interaction_list

    def _process_frame_wrapper(self, frame_idx):
        """
        Wrapper for the MD Trajectory procession.

        Parameters
        ----------
        frame_idx : int
            Number of the frame to be processed.

        Returns
        -------
        tuple
            Tuple containing the frame index and the result of from the process_frame function.
        """
        return frame_idx, self._process_frame(frame_idx)

    def _fill_missing_frames(self, df):
        """
        Fills the frames with no interactions in the DataFrame with placeholder values.

        Parameters
        ----------
        df : pd.DataFrame 
            The input DataFrame with frames that have no interactions.

        Returns
        -------
        pd.DataFrame 
            DataFrame with placeholder values in the frames with no interactions.
        """

        # Create a set containing all unique values in the 'FRAME' column
        existing_frames = set(df["FRAME"])

        # Create a list to store new rows for missing numbers
        missing_rows = []

        # Iterate through numbers from 0 to md_len
        for frame_number in range(1, self.md_len):
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

    def _process_trajectory(self):
        """
        Process protein-ligand trajectory with multiple CPUs in parallel.

        Returns
        -------
        pd.DataFrame 
            A DataFrame containing all the protein-ligand interaction data from the whole trajectory.
        """
        if self.dataframe is None:
            print("\033[1mProcessing protein-ligand trajectory\033[0m")
            print(f"\033[1mUsing {self.num_processes} CPUs\033[0m")

            with Pool(processes=self.num_processes) as pool:
                frame_indices = range(1, self.md_len)

                # Initialize the progress bar with the total number of frames
                pbar = tqdm(
                    total=self.md_len - 1,
                    ascii=True,
                    desc="\033[1mAnalyzing frames\033[0m",
                )

                results = []
                for result in pool.imap(self._process_frame_wrapper, frame_indices):
                    results.append(result)
                    pbar.update(1)  # Update the progress manually

            # Close the progress bar
            pbar.close()

            # Extract the results and sort them by frame index
            results.sort(key=lambda x: x[0])
            interaction_lists = [result[1] for result in results]

            interaction_list = pd.concat(interaction_lists)

            interaction_list.to_csv("interactions_gathered.csv")

        elif self.dataframe is not None:
            print(f"\033[1mGathering data from {self.dataframe}\033[0m")
            interaction_tmp = pd.read_csv(self.dataframe)
            interaction_list = interaction_tmp.drop(interaction_tmp.columns[0], axis=1)

        interaction_list["Prot_partner"] = (
            interaction_list["RESNR"].astype(str)
            + interaction_list["RESTYPE"]
            + interaction_list["RESCHAIN"]
        )

        interaction_list = self._fill_missing_frames(
            interaction_list,
        )

        print("\033[1mProtein-ligand trajectory processed\033[0m")

        return interaction_list
