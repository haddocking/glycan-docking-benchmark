#!/usr/bin/env python3
#PBS -N glycan_res
#PBS -q medium
#PBS -l nodes=1:ppn=48
#PBS -S /bin/python

# Add the glycan residues to the info text files
# Execution: python add_glycan_res.py location/of/the/pdb/files
import argparse
import os

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

    for pdbf in files_list:
        if pdbf.endswith("_ref.pdb"):
            pass
        elif pdbf.endswith("_r_b.pdb"):
            pass
        elif pdbf.endswith("_r_u.pdb"):
            pass
        elif pdbf.endswith("_info.txt"):
            pass
        elif pdbf.endswith("_l_b.pdb"):
            ligand_file = open(pdbf, "rt")
            residues_list = []
            for line in ligand_file:
                splitted_line = line.split()
                id = splitted_line[0]
                if id == "HETATM":
                    residues = splitted_line[3]
                    residues_list.append(residues)
                    residues_list = sorted(set(residues_list))
            residues_found = map(str, residues_list)
            residues_found = ', '.join(residues_found)
            ligand_file.close()

    for txtf in files_list:
        if txtf.endswith("_ref.pdb"):
            pass
        elif txtf.endswith("_r_b.pdb"):
            pass
        elif txtf.endswith("_l_b.pdb"):
            pass
        elif txtf.endswith("_info.txt"):
            info_file = open(txtf, "a")
            info_file.write(residues_found + os.linesep)
            info_file.close()
