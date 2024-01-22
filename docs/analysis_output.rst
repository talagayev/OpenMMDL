**OpenMMDL Analysis Output Files**
=================================

The **OpenMMDL Analysis** script creates output folders and files during and after the analysis.

The following list contains an detailed overview of the folders and files:
The Optional files are highlighted with an *asteriks**

Barcodes
------------------------------
**Barcodes**: A folder that contains the figures of the barcodes for the interactions.



.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Name
     - Description
   * - acceptor_interactions.png*
     - File displaying the barcodes for all hydrogen bond acceptor interactions.
   * - donor_interactions.png*
     - File displaying the barcodes for all hydrogen bond donor interactions.
   * - halogen_interactions.png*
     - File displaying the barcodes for all halogen interactions.
   * - hydrophobic_interactions.png*
     - File displaying the barcodes for all hydrophobic interactions.
   * - metal_interactions.png*
     - File displaying the barcodes for all metal interactions.
   * - saltbridge_interactions.png*
     - File displaying the barcodes for all saltbridge interactions.
   * - pication_interactions.png*
     - File displaying the barcodes for all pication interactions.
   * - pistacking_interactions.png*
     - File displaying the barcodes for all pistacking interactions.
   * - waterbridge_interactions.png*
     - File displaying the barcodes for all waterbridge interactions.

Binding_Modes_Markov_States
------------------------------

**Binding_Modes_Markov_States**: A folder that contains the figures of the Markov state figures and 2D depiction figure of the  top 10 occuring binding modes.


.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Name
     - Description
   * - checkpoint.chk
     - Checkpoint saved every 10000 steps.
   * - 10x_checkpoint.chk
     - Checkpoint saved every 100000 steps.
   * - 100x_checkpoint.chk*
     - Checkpoint saved every 1000000 steps.

RMSD
------------------------------
**RMSD**: A folder that contains the files for the RMSD calculation during the simulation. The RMSD over time is default, while the RMSD between the frames is optional

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Name
     - Description
   * - RMSD_over_time.csv
     - CSV File with the values of the RMSD  of the protein/backbone/ligand during the simulation time.
   * - RMSD_over_time.png
     - Figure of the RMSD of the protein/backbone/ligand during the simulation time.
   * - RMSD_between_the_frames.png*
     - Matrix figure displaying the RMSD between each consecutive frame.


*Minimization Equilibration*: Folder that contains topology files after the minimization and equilibration with **OpenMM**.

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Name
     - Description
   * - Energyminimization_protein_name.pdb
     - Prepared PDB File after OpenMM energy minimization.
   * - Equilibration_protein_name.pdb
     - Prepared PDB File after OpenMM energy minimization and equilibration.

*MD Output*: Folder that contains the Output trajectory files generated during the MD Simulation.


.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Name
     - Description
   * - output_protein_name.pdb
     - PDB Topology File of the first frame of the simulation.
   * - trajectory.dcd
     - Trajectory of the OpenMM Simulation.


MD Postprocessing
------------------------------
**MD Postprocessing**: A folder that contains the postprocessing files after the MD Simulation.

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Name
     - Description
   * - centered_old_coordinates_top.pdb
     - Topology File of the centered protein in PDB Format.
   * - centered_old_coordinates_top.gro
     - Topology File of the centered protein in Gromacs GRO Format.
   * - centered_old_coordinates.dcd
     - Trajectory File of the centered protein in DCD Format.
   * - centered_old_coordinates.xtc
     - Trajectory File of the centered protein in XTC Format.
   * - centered_traj_unaligned.dcd
     - Trajectory File of the unaligned centered protein in DCD Format with all atoms and new coordinates.
   * - centered_traj_unaligned.xtc
     - Trajectory File of the unaligned centered protein in XTC Format with all atoms and new coordinates.
   * - prot_lig_traj_unaligned.dcd*
     - Trajectory File of the unaligned centered protein in DCD Format with only protein and ligand atoms and new coordinates.
   * - prot_lig_traj_unaligned.xtc*
     - Trajectory File of the unaligned centered protein in XTC Format with only protein and ligand atoms and new coordinates.
     
Final Output
------------------------------
**Final Output**: A folder that contains the final files after the MD Simulation, ready to be analyzed. This folder contains 2 subfolders.

1. *All Atoms*: Folder that contains the centered topology files of all atoms with new coordinates according to the center of mass.

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Name
     - Description
   * - centered_top.pdb*
     - Topology File of the centered protein in PDB Format with all atoms and new coordinates.
   * - centered_top.gro*
     - Topology File of the centered protein in Gromacs GRO Format with all atoms and new coordinates.
   * - centered_traj.dcd*
     - Trajectory File of the aligned centered protein in DCD Format with all atoms and new coordinates.
   * - centered_traj.xtc*
     - Trajectory File of the aligned centered protein in XTC Format with all atoms and new coordinates.



2. *Prot Lig*: Folder that contains the centered topology files of only the protein and ligand atoms with new coordinates according to the center of mass.



.. list-table::
   :header-rows: 1
   :widths: 25 75


   * - Name
     - Description
   * - prot_lig_top.pdb*
     - Topology File of the centered protein in PDB Format with only protein and ligand atoms and new coordinates.
   * - prot_lig_top.gro*
     - Topology File of the centered protein in Gromacs GRO Format with only protein and ligand atoms and new coordinates.
   * - prot_lig_traj.dcd*
     - Trajectory File of the centered protein in DCD Format with only protein and ligand atoms and new coordinates.
   * - prot_lig_traj.xtc*
     - Trajectory File of the centered protein in XTC Format with only protein and ligand atoms and new coordinates.
