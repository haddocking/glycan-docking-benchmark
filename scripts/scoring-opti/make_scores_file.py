#!/usr/bin/python

import argparse
import sys
import os
import re
import numpy as np
import csv

'''Given the PDBs and the i-RMSD.dat, this script generates the .scores file'''

#complex = sys.argv[1]
#path = complex + 'run-true-interface/structures/it0/'

weights = [1.0, 1.0, 1.0, 1.0, 1.0]


def get_haddock_terms(file_nam, ranking):

    with open(ranking, 'r') as inp_fh: 
            rmsd = {str('{}_{}'.format(lines.split('\t')[0], lines.split('\t')[1])): float(lines.split('\t')[2]) for lines in inp_fh} 
            
    with open(file_nam, 'r') as inp_fh:
        os.chdir(path)
        haddock_terms = []
        for line in inp_fh:
            if len(line) > 4:
                rms = rmsd[line.split('.')[0]]
                structure = line.strip().split()[0]               
                with open(structure, 'r') as pdb_fh:
                    for line in pdb_fh:
                        if line[0:6] == 'REMARK':
                            if re.match(r'REMARK energies', line):
                                ene = line.split(',')
                                Eelec = ene[6].strip()
                                Evdw = ene[5].strip()
                                Eair = ene[7].strip()
                            elif re.match(r'REMARK Desolvation energy', line):
                                Edesolv = line.split()[3].strip()
                            elif re.match(r'REMARK buried surface area', line):
                                BSA = line.split()[4].strip()
                            else:
                                continue
                    energy = [float(component) for component in [Eelec, Evdw, Edesolv, Eair, BSA]]
                    hs = np.dot(energy, weights)
                    terms = [structure, rms, str(hs), Eelec, Evdw, Edesolv, Eair, BSA]
                    haddock_terms.append(terms)
    return haddock_terms


parser = argparse.ArgumentParser()
parser.add_argument("pdbfiles_folders", help="PDB files folders location")
args = parser.parse_args()
os.chdir(args.pdbfiles_folders)

curr_dir = os.getcwd()
directories_list = []
for direct in os.listdir():
    if os.path.isdir(direct):
        directories_list.append(direct)

for pdb_id in directories_list:
    w_direct = curr_dir + '/' + pdb_id + '/'
    os.chdir(w_direct)
    path = 'run-true-interface/structures/it0/'

    if __name__ == '__main__':

        print('Extracting haddock terms')

        haddock_terms = get_haddock_terms(path + 'file.nam', 'i-RMSD.dat')
        #print(haddock_terms)
        matrix = sorted(haddock_terms, key=lambda tup: tup[2], reverse=True)
        #print(matrix)

        print('Creating file')

        with open('/trinity/login/apeliss/glycan-docking/05-score-opti/true-interface/' + pdb_id + '.scores', 'w') as out_fh:
            writer = csv.writer(out_fh, quoting=csv.QUOTE_MINIMAL)
            for row in matrix:
                writer.writerow(row)
