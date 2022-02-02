#!/usr/bin/env python3
#PBS -N dataset
#PBS -q medium
#PBS -l nodes=1:ppn=48
#PBS -S /bin/python

# Create folders, for each complex, containing the original complex,
# receptor and glycan (ligand), bound & unbound structures pdb files
# and an info.txt file containing the glycan name and residues.

# Execution: python organize_benchmark.py location/of/the/pdb/files
import pandas as pd
import argparse
import shutil
import os
import re

parser = argparse.ArgumentParser()
parser.add_argument("dataset_folder", help="Dataset folder location")
args = parser.parse_args()
os.chdir(args.dataset_folder)
target_list = os.listdir()

dataset = pd.read_excel('gd-dataset.xlsx')
dataset.to_csv('dataset.csv', index=None, header=True)
pd.set_option('max_colwidth', 100)

# Create the name of the folders for each complex: PDBID
folders = []
for structure in target_list:
    pattern = re.compile(r'_ref\.pdb$')
    matches = pattern.findall(structure)
    if matches:
        for match in matches:
            selected = structure[0:4]
            folders.append(selected)

# Create the information text files (PDBID_info.txt) containing the glycan name
for pdb_id in folders:
    glycan_name = str(dataset[dataset['PDB_ID'] == pdb_id]['Glycoligand_Sequence'])
    glycan_name = glycan_name[5:-42]
    with open(pdb_id + '_info.txt', 'w') as info_file:
        info_file.write('Glycan name: ' + glycan_name)
        info_file.write('\nResidues: ')

# Clean the receptor unbounded structure (PDBID_r_u.pdb),
# and create the receptor bounded structure (PDBID_r_b.pdb)
# and ligand bounded structure pdb files (PDBID_l_b.pdb)
for pdbf in target_list:
    if pdbf.endswith('_r_u_ori.pdb'):
        with open(pdbf[0:8] + '.pdb', 'w') as unbound:
            with open(pdbf, 'rt') as unbound_ori:
                pattern = re.compile('^ATOM|^TER|^END')
                for line in unbound_ori:
                    if re.search(pattern, line):
                        unbound.write(line)
        os.remove(pdbf)

structures_ref = []
for pdbf in target_list:
    if pdbf.endswith('_ref.pdb'):
        structures_ref.append(pdbf)

for ref in structures_ref:
    with open(ref[0:4] + '_r_b.pdb', 'w') as receptor:
        with open(ref, 'rt') as ref:
            pattern = re.compile('^ATOM|^TER|^END')
            for line in ref:
                if re.search(pattern, line):
                    receptor.write(line)

for ref in structures_ref:
    with open(ref[0:4] + '_l_b.pdb', 'w') as ligand:
        with open(ref, 'rt') as ref:
            pattern = re.compile('^HETATM|^END')
            for line in ref:
                if re.search(pattern, line):
                    ligand.write(line)

# Create folders for each PDBID structure and
# move each file to its corresponding PDBID folder
all_files = os.listdir()

for pdb_id in folders:
    if not os.path.isdir(pdb_id):
        os.mkdir(pdb_id)

for structure in all_files:
    for pdb_id in folders:
        if structure[0:4] == pdb_id:
            shutil.move(structure, pdb_id)

for pdb_id in folders:
    unbound = str(dataset[dataset['PDB_ID'] == pdb_id]['Unbound_PDB'])
    unbound = unbound[5:-33].strip()
    for structure in all_files:
        if structure.endswith('_r_u.pdb'):
            if structure[0:4] == unbound:
                shutil.copy(structure, pdb_id)

for unbound in all_files:
    if unbound.endswith('_r_u.pdb'):
        os.remove(unbound)
