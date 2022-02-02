## Glycan-docking-benchmark

The repository contains the following information:

# HADDOCK-ready

Directory composed by sub-directories corresponding to each complex of the glycan benchmark. 


Each sub-directory contains HADDOCK-ready files:

* 'XXXX_r_b.pdb' : Bound receptor PDB
* 'XXXX_l_b.pdb' : Bound ligand PDB
* 'XXXX_r_u.pdb' : Unbound receptor PDB
* 'XXXX_l_u.pdb' : Unbound ligand PDB

The original file:
* 'XXXX_ref.pdb' : Original complex structure PDB

An info file:
* 'XXXX_info.txt' : Associated information file for each PDB complex containing the glycan name and its residue names

And the distance restraints file:

* 'true_interface.tbl' : Ambiguous interactions restraints based on true-interface restraints


Each sub-directory contains an analysis file containing the following files:

* 'XXXX_analysis.pdb' : Processed reference complex PDB file to make it ready for analysis with ProFit (v3.3)
