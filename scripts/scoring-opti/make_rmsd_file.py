# From the rmsd.csv files, creates an i-RMSD.dat file for each complex
# which is a table with rmsd data of all the PDB docked structures from the docking:
# fieldnames --> pdb_id, rank, rmsd
# Usage: python3 make_rmsd_file.py location_of_pdb_files/

import pandas as pd
import argparse
import glob
import csv
import os

parser = argparse.ArgumentParser()
parser.add_argument("pdbfiles_folders", help="PDB files folers location")
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
    with open('rmsd.csv', 'r') as csvfile:
        # rmsd_reader = csv.reader(csvfile, delimiter=',')
        headers = ['pdb_id', 'iteration', 'rank',
                   'i-rmsd', 'fnat', 'l-rmsd', 'i-l-rmsd',
                   'haddock_score', 'path_to_the_structure']
        rmsd_csv = pd.read_csv(csvfile, names=headers)
        rmsd_csv = rmsd_csv[rmsd_csv['iteration'] == 'it0']
        rmsd_data = rmsd_csv[['pdb_id', 'rank', 'i-rmsd']]
        rmsd_data_no_headers = rmsd_data.rename(columns=rmsd_data.iloc[0]).drop(rmsd_data.index[0])
        rmsd_data_no_headers.to_csv('i-RMSD.dat', index=False, sep='\t')
        # rmsd_data_string = rmsd_data_no_headers.to_string(index=False)
        # with open('i-RMSD.dat', 'w', newline='') as datfile:
        #     dat_file_csv = csv.writer(datfile, delimiter='\t')
        #     for row in rmsd_data_no_headers:
        #         dat_file_csv.writerow(row)
