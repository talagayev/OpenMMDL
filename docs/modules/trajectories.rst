API Documentation for TrajectorySaver
=========================================

.. py:class:: TrajectorySaver(pdb_md, ligname, special, nucleic)

    Saves frames and subsets of molecular dynamics trajectories involving ligands, receptors, and interacting waters.

    :param mda.Universe pdb_md: MDAnalysis Universe object containing the trajectory and topology.
    :param str ligname: Name of the ligand in the PDB topology.
    :param str special: Name of any special residue (e.g. HEM) to include.
    :param bool nucleic: Whether the receptor is nucleic acid (`True`) or protein (`False`).


    .. py:method:: save_interacting_waters_trajectory(interacting_waters, outputpath)

        Saves `.pdb` and `.dcd` trajectory files including the ligand, receptor, and specified interacting water molecules.

        :param list interacting_waters: List of water residue IDs to include in the output.
        :param str outputpath: Path where the output `.pdb` and `.dcd` files will be written. Defaults to `'./Visualization/'`.

        :returns: Writes files `interacting_waters.pdb` and `interacting_waters.dcd` to the specified path.
        :rtype: None


    .. py:method:: save_frame(frame, outpath, selection=False)

        Saves a single frame from the trajectory to file, with optional atom selection.

        :param int frame: Index of the frame to save.
        :param str outpath: File path to write the frame to.
        :param str selection: Optional MDAnalysis selection string. If not provided, all atoms are saved. Defaults to `False`.

        :returns: Saves the selected frame as a structure file (e.g., `.pdb`).
        :rtype: None
