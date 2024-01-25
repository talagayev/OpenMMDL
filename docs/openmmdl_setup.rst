**Running OpenMMDL Setup**
=============================

This page displays the preparation paths and forcefields available in **OpenMMDL** and showcases the application of **OpenMMDL Setup**.

.. figure:: /_static/images/OpenMMDL_Setup.png
    :figwidth: 600px
    :height: 100px
    :align: center
|  
To start the **OpenMMDL Setup** we need to activate the `openmmdl` environment. To do this we have to enter the following command lines:

.. code-block:: text

    conda activate openmmdl

Now that we have activated the `openmmdl` environment we can start **OpenMMDL Setup**. To do this you need to type the following:

.. code-block:: text

    openmmdl_setup

This will open the **OpenMMDL Setup**, which you can use for the creation of the input files for **OpenMMDL Simulation**.

There are two possible options to create the input files for **OpenMMDL Simulation**:

1. The PDBFixer path, where a `pdb` file of the protein is used as an input for the preparation and simulation.
The tutorial for the PDBFixer path can be found :doc:`here </tutorial_pdb_path>`.

Here is the table of the currently available forcefields and water models for the PDBFixer path: 

.. figure:: /_static/images/Forcefield_watermodels.png
   :figwidth: 725px
   :align: center

2. The Amber path, where `prmtop` and `inpcrd` files are used for the preparation and simulation. This path allows us to either use already prepared `prmtop` and `inpcrd` as an input or create the `prmtop` and `inpcrd` from PDB files of the receptor and ligand.
The tutorial for the Amber path can be found :doc:`here </tutorial_amber_path>`.

.. figure:: /_static/images/amber_ff.png
   :figwidth: 725px
   :align: center

In the table, the first row is the default setting, and the term `other` allows users to type their desired forcefields from those accessible in AmberTools 22.0 into the designated textbox.
