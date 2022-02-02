# Importing the script calc_distances.py, calculate the pairwise distances between all atoms of the 
# receptor and the ligand of a complex pdb file of the receptor and ligand structures pdb files 
# found in their folder named by their PDBID from the original structure file.

# Execution: python iteration_calc-distances.py location/of/the/folders-for-each-structure-with-name-PDBID, e.g., python .\iteration_calc-distances.py .\HADDOCK-ready\ 
import os
import argparse
import calc_distances

parser = argparse.ArgumentParser()
parser.add_argument("pdbfiles_folders", help="PDB files folders location")
args = parser.parse_args()
os.chdir(args.pdbfiles_folders)

curr_dir = os.getcwd()
directories_list = []
for subdir, dirs, files in os.walk(curr_dir):
    for direct in dirs:
        directories_list.append(direct)

for structure_direct in directories_list:
    w_direct = curr_dir + '/' + structure_direct + '/'
    os.chdir(w_direct)

 
    files_list = []
    
    for filename in os.listdir():
        files_list.append(filename)

    
    for receptor_file in files_list:
        if receptor_file.endswith("r_b.pdb"):
            for ligand_file in files_list:
                if ligand_file.endswith("l_b.pdb"):
                    calc_distances.interface_residues(receptor_file, ligand_file)
