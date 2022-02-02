# This script creates a table stating the following:
# fieldnames -> PDB ID complex, protein classification, experimental method,
#               resolution, number of glycan sugar units, linearity
# Use in folder where all the pdb information files and the pdb structure
# files are found.
# Usage: python3 scripts/table-complete-gd.py 01-organization/05-gd-dataset

import pandas as pd
import collections
import argparse
import os
import re

parser = argparse.ArgumentParser()
parser.add_argument("pdbfiles", help="PDB files folder location")
args = parser.parse_args()
os.chdir(args.pdbfiles)

# List with reference complex pdb files
ref_complexes = []
for complex in os.listdir():
    if complex.endswith('_ref.pdb'):
        ref_complexes.append(complex)

# List with glycan-bound files
glycan_list = []
for glycan in os.listdir():
    if glycan.endswith('_l_b.pdb'):
        glycan_list.append(glycan)

# Create a list with all PDB IDs of each structure
pdb_IDs = []
for complex in ref_complexes:
    pdb_ID = complex[0:4]
    pdb_IDs.append(pdb_ID)
pdb_IDs = sorted(pdb_IDs)

# Create a dictionary where keys are the PDB ID of each structure, and values
# are its protein classification.
prot_class = {}
for complex in ref_complexes:
    with open(complex, 'r') as fh:
        for line in fh.readlines():
            if line.startswith('HEADER'):
                prot = str(line[10:50])
                prot = prot.rstrip()
                for name in pdb_IDs:
                    if name == complex[0:4]:
                        prot_class[name] = prot
prot_class = collections.OrderedDict(sorted(prot_class.items()))

# Create a dictionary where keys are the PDB ID of each structure, and values
# are the number of units of the sugar the protein is interacting with.
sugar_units = {}
for glycan in glycan_list:
    with open(glycan, 'r') as fh:
        for line in fh.readlines():
            if line.startswith('HETATM'):
                residue = int(line[23:26])
        for name in pdb_IDs:
            if name == glycan[0:4]:
                sugar_units[name] = residue
sugar_units = collections.OrderedDict(sorted(sugar_units.items()))

# Create a dictionary where keys are the PDB ID of each structure, and values
# are the experimental method used to determine their structure.
exp_method = {}
for complex in ref_complexes:
    with open(complex, 'r') as fh:
        for line in fh.readlines():
            if line.startswith('EXPDTA'):
                prot = str(line[10:50])
                prot = prot.rstrip()
                for name in pdb_IDs:
                    if name[0:4] == complex[0:4]:
                        exp_method[name] = prot
exp_method = collections.OrderedDict(sorted(exp_method.items()))

# Create a dictionary where keys are the PDB ID of each structure, and values
# are the resolution of the experimental structure.
resolution = {}
for complex in ref_complexes:
    with open(complex, 'r') as fh:
        for line in fh.readlines():
            if line.startswith('REMARK   2 RESOLUTION.'):
                prot = str(line[26:30])
                prot = prot.rstrip()
                for name in pdb_IDs:
                    if name[0:4] == complex[0:4]:
                        resolution[name] = prot
resolution = collections.OrderedDict(sorted(resolution.items()))

# Create a dictionary where keys are the PDB ID of each structure, and values
# are the linearity of each structure's saccharide.
dataset = pd.read_excel('gd-dataset.xlsx')
dataset.to_csv('dataset.csv', index=None, header=True)
pd.set_option('max_colwidth', 100)

glycans_glycam = {}
for pdbid in dataset['PDB_ID']:
    sequence = str(dataset[dataset['PDB_ID'] == pdbid]['Glycoligand_Sequence'])
    sequence = sequence[5:].strip()
    sequence = sequence.partition('\n')
    glycans_glycam[pdbid] = sequence[0]

linearity = {}
branched = []
linear = []

for pdbid, sequence in glycans_glycam.items():
    pattern = re.compile(r'\[')
    matches = pattern.findall(sequence)
    if matches:
        for name in pdb_IDs:
            if name == pdbid:
                linearity[name] = 'branched'
    else:
        for name in pdb_IDs:
            if name == pdbid:
                linearity[name] = 'linear'
linearity = collections.OrderedDict(sorted(linearity.items()))


with open('table', 'w') as tbl:
    tbl.write('PDB_ID' + '\t' + 'Protein_class' + '\t' +
              'Experimental_method' + '\t' + 'Resolution' +
              '\t' + 'Linearity' + '\t' + 'Number_of_sugar_units' + '\n')
    for (i, j, k, l, m, n) in zip(pdb_IDs, prot_class, exp_method,
                                  resolution, linearity, sugar_units):
        tbl.write(i + '\t' + prot_class[j] + '\t' + exp_method[k] + '\t' +
                  str(resolution[l]) + '\t' + linearity[m] + '\t' +
                  str(sugar_units[n]) + '\n')
