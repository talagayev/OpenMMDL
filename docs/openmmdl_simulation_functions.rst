OpenMMDL simulation functions
=============================

This page displays all the functions of **OpenMMDL Simulation**.

openmmdl_simulation.scripts.cleaning_procedures
------------------------------

.. py:function:: cleanup(protein_name)
    
    Cleans up the PDB Reporter Output File and MDTraj Files of the performed simulation.
    
    :param str protein_name: Name of the protein PDB.

    :returns: None.
    :rtype: None
   

.. py:function:: create_directory_if_not_exists(directory_path)
    
    Create a directory if it doesn't exist, or overwrite it if already does.
    
    :param str directory_path: Path of the directory that you want to create.

    :returns: None.
    :rtype: None


.. py:function:: copy_file(src, dest)
    
    Copy a file to the destination path.
    
    :param str src: Path of the file that needs to be copied.
    :param str dest: Path of destination where the file needs to be copied to.

    :returns: None.
    :rtype: None


.. py:function:: organize_files(source, destination)
    
    Organizes the files and moves them from the source to the destination directory.
    
    :param str source: Path of the file that needs to be moved.
    :param str destination: Path of destination where the file needs to be moved to.

    :returns: None.
    :rtype: None




.. py:function:: post_md_file_movement(protein_name: str, prmtop: str = None, inpcrd: str = None, ligands: List[str] = None)
    
    Organizes and moves the files after the MD simulation to their respective directories.
    
    :param str protein_name: Name of the protein PDB.
    :param prmtop: Path to the AMBER topology file.
    :param inpcrd: Path to the AMBER coordinate file.
    :param ligands: List of paths to the ligand files.
    :type prmtop: Optional [str]
    :type inpcrd: Optional [str]
    :type ligands: Optional [List[str]]


    :returns: None.
    :rtype: None
