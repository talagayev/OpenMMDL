import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

import MDAnalysis as mda
import MDAnalysis.transformations as trans

from MDAnalysis.analysis import rms, diffusionmap, align
from MDAnalysis.analysis.distances import dist

def rmsd_for_atomgroups(prot_lig_top_file, prot_lig_traj_file, selection1, selection2=None):
    """Calulate the RMSD for selected atom groups, and save the csv file and plot.

    Parameters
    ----------
    prot_lig_top_file: str
        name of the PDB input file.
    prot_lig_traj_file: str
        name of the DCD input file.
    selection1: str
        Selection string for main atom group, also used during alignment.
    selection2: list of str, optional
        Selection strings for additional atom groups.

    Returns
    -------
    rmsd_df: pandas.core.frame.DataFrame
        DataFrame containing RMSD of the selected atom groups over time.
    
    """ 
    universe=mda.Universe(prot_lig_top_file, prot_lig_traj_file)    
    universe.trajectory[0]
    ref = universe
    rmsd_analysis = rms.RMSD(universe, ref, select=selection1, groupselections=selection2)
    rmsd_analysis.run()
    columns = [selection1, *selection2] if selection2 else [selection1]
    rmsd_df = pd.DataFrame(np.round(rmsd_analysis.rmsd[:, 2:], 2), columns=columns)
    rmsd_df.index.name = "frame"

    rmsd_df.to_csv('./RMSD/RMSD_over_time.csv', sep=' ')

    rmsd_df.plot(title="RMSD of protein and ligand")
    plt.ylabel("RMSD (Å)")
    plt.savefig('./RMSD/RMSD_over_time.png')

    return rmsd_df

def RMSD_dist_frames(prot_lig_top_file, prot_lig_traj_file, lig, nucleic=False):
    """Calculate the RMSD between all frames in a matrix.

    Parameters
    ----------
    prot_lig_top_file: str
        name of the PDB file.
    prot_lig_traj_file: str
        name of the DCD file.
    lig: str
        ligand name saved in the above pdb file. Selection string for the atomgroup to be investigated, also used during alignment.

    Returns
    -------
    pairwise_rmsd_prot: np.ndarray
        Numpy array of RMSD values for pairwise protein structures.
    pairwise_rmsd_lig: np.ndarray
        Numpy array of RMSD values for ligand structures.

    """
    universe=mda.Universe(prot_lig_top_file, prot_lig_traj_file)
    if nucleic:
        pairwise_rmsd_prot = diffusionmap.DistanceMatrix(universe, select="nucleic").run().dist_matrix
    else:
        pairwise_rmsd_prot = diffusionmap.DistanceMatrix(universe, select="protein").run().dist_matrix
    pairwise_rmsd_lig = diffusionmap.DistanceMatrix(universe, f"resname {lig}").run().dist_matrix

    max_dist = max(np.amax(pairwise_rmsd_lig), np.amax(pairwise_rmsd_prot))
    
    fig, ax = plt.subplots(1,2)
    fig.suptitle("RMSD between the frames")

    # protein image
    img1 = ax[0].imshow(pairwise_rmsd_prot, cmap="viridis", vmin=0, vmax=max_dist)
    if nucleic:
        ax[0].title.set_text("nucleic")
    else:
        ax[0].title.set_text("protein")
    ax[0].set_xlabel("frames")
    ax[0].set_ylabel("frames")
    
    # ligand image
    img2 = ax[1].imshow(pairwise_rmsd_lig, cmap="viridis", vmin=0, vmax=max_dist)
    ax[1].title.set_text("ligand")
    ax[1].set_xlabel("frames")

    fig.colorbar(img1, ax=ax, orientation="horizontal", fraction=0.1, label="RMSD (Å)")

    plt.savefig('./RMSD/RMSD_between_the_frames.png')
    return pairwise_rmsd_prot, pairwise_rmsd_lig
