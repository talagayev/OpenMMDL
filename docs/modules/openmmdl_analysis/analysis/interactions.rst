API Documentation for interactions
============================


.. py:class:: InteractionAnalyzer(pdb_md, dataframe, num_processes, lig_name, special_ligand, peptide, md_len)

    Analyzes molecular interactions between a protein and a ligand/peptide
    throughout an MD trajectory using ProLIF.

    :ivar mda.Universe pdb_md: MDAnalysis Universe object representing the topology and trajectory.
    :ivar str or None dataframe: Path to an existing interaction CSV file. If None, the trajectory will be processed anew.
    :ivar int num_processes: Number of CPU cores to use for parallel frame analysis.
    :ivar str lig_name: Residue name of the ligand in the complex.
    :ivar str special_ligand: Residue name for special ligands like metal ions (optional).
    :ivar str peptide: Chain ID of the peptide ligand (optional).
    :ivar int md_len: Number of frames in the trajectory.
    :ivar pd.DataFrame interaction_list: DataFrame storing the extracted interactions across the trajectory.

    .. py:method:: _select_ligand_group(self, frame)

        Determine the MDAnalysis atom selections corresponding to the ligand and
        receptor for a given frame.

        :param int frame: Frame index to extract atom selections from.
        :returns: Tuple containing the ligand and receptor atom groups.
        :rtype: tuple

    .. py:method:: _run_prolif(self, frame)

        Executes the ProLIF fingerprint calculation for the requested frame.

        :param int frame: Frame index to evaluate.
        :returns: Raw ProLIF dataframe containing interaction information.
        :rtype: pd.DataFrame

    .. py:method:: _convert_prolif_dataframe(self, prolif_df, frame)

        Convert the raw ProLIF dataframe into the legacy PLIP-like schema used by
        the rest of the analysis pipeline.

        :param pd.DataFrame prolif_df: Raw results produced by ProLIF.
        :param int frame: Frame index associated with the dataframe.
        :returns: DataFrame with information retrieved from ProLIF.
        :rtype: pd.DataFrame

    .. py:method:: _process_frame(self, frame)

        Process a single frame of MD simulation.

        :param int frame: The number of the frame that will be processed.
        :returns: A dataframe conatining the interaction data for the processed frame.
        :rtype: pd.DataFrame

    .. py:method:: _process_frame_wrapper(self, frame_idx)

        Wrapper for the MD Trajectory procession.

        :param int frame_idx: Number of the frame to be processed.
        :returns: Tuple containing the frame index and the result of from the process_frame function.
        :rtype: tuple

    .. py:method:: :fill_missing_frames(self, df)

        Fills the frames with no interactions in the DataFrame with placeholder values.

        :param pd.DataFrame df: The input DataFrame with frames that have no interactions
        :returns: DataFrame with placeholder values in the frames with no interactions.
        :rtype: pd.DataFrame

    .. py:method:: _process_trajectory(self)

        Process protein-ligand trajectory with multiple CPUs in parallel.

        :returns: A DataFrame containing all the protein-ligand interaction data from the whole trajectory.
        :rtype: pd.DataFrame
