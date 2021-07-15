## Glycan-docking-benchmark

The repository contains the following information:

# HADDOCK-ready-files

Directory composed by sub-directories corresponding to each complex of the glycan benchmark. 

Each sub-directory contains HADDOCK-ready files:

* 'XXXX_ref.pdb' : Original complex structure PDB
* 'XXXX_r_b_ori.pdb' :
* 'XXXX_l_b_ori.pdb' :
* 'XXXX_r_b.pdb' :
* 'XXXX_l_b.pdb' :

* 'XXXX_info.txt' : Associated information file for each PDB complex containing the glycan name and its residue names.

The files containing a list of the interface residues extracted from the bound complex (residues whose atoms are within a 3.9 A cut-off of any atom of any partner):

* 'ligand_con.txt' : Ligand interface residues list
* 'receptor_con.txt' : Receptor interface residues list

And the distance restraints file:

* 'true_interface.tbl' : Ambiguous interactions restraints based on true-interface restraints.
