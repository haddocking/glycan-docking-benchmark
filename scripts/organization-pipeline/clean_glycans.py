#!/usr/bin/env python3
#PBS -N clean_glycans
#PBS -q medium
#PBS -l nodes=1:ppn=48
#PBS -S /bin/python

# Removes the residues in the glycans pdb files (PDBID_l_b.pdb) that are not
# belonging to sugars (e.g., HOH).

# Execution: python clean_glycans.py location/of/the/pdbid/folders
import argparse
import glob
import os
import re

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
            with open(pdbf[0:-4], 'w') as new_pdbf:
                with open(pdbf, 'rt') as pdbf:
                    pattern1 = re.compile('^HETATM')
                    pattern2 = re.compile('^END')
                    for line in pdbf:
                        if re.search(pattern1, line):
                            resnum = int(line[22:26])
                            if resnum < 12:            # max resnum is 11 in this specific dataset
                                new_pdbf.write(line)
                        elif re.search(pattern2, line):
                            new_pdbf.write(line)

    pre_files = glob.glob('*_l_b.pdb')
    post_files = glob.glob('*_l_b')
    for file in pre_files:
        os.unlink(file)
    for file in post_files:
        os.rename(file, file + '.pdb')
