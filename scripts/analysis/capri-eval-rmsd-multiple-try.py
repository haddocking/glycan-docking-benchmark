#!/usr/bin/env python3

import subprocess
import argparse
import glob
import csv
import os


def get_rank_and_hs(file_list):
    # file_list = 'run-true-interface/structures/it0/file.list'
    output_dict = {}
    with open(file_list, 'r') as fh:
        for i, line in enumerate(fh.readlines(), start=1):
            haddock_score = float(line.split()[-2])
            pdb_name = line.split()[0].split(':')[-1][:-1]
            rank = i
            output_dict[pdb_name] = {'ranking': rank,
                                     'haddock-score': haddock_score}
    return output_dict


parser = argparse.ArgumentParser()
parser.add_argument("pdbfiles_folders", help="PDB files folders location")
args = parser.parse_args()
os.chdir(args.pdbfiles_folders)

curr_dir = os.getcwd()
directories_list = []
for direct in os.listdir():
    directories_list.append(direct)

for pdb_id in directories_list:
    w_direct = curr_dir + '/' + pdb_id + '/'
    os.chdir(w_direct)
    ref_pdb = pdb_id + '_profit-modif.pdb'
    pdb_it0 = glob.glob('run-true-interface/structures/it0/*pdb')
    pdb_it1 = glob.glob('run-true-interface/structures/it1/*pdb')
    pdb_itw = glob.glob('run-true-interface/structures/it1/water/*w.pdb')

    file_list_it0 = 'run-true-interface/structures/it0/file.list'
    file_list_it1 = 'run-true-interface/structures/it1/file.list'
    file_list_itw = 'run-true-interface/structures/it1/water/file.list'

    score_ranking_it0 = get_rank_and_hs(file_list_it0)
    score_ranking_it1 = get_rank_and_hs(file_list_it1)
    score_ranking_itw = get_rank_and_hs(file_list_itw)

    rmsd_it0_data = []
    for path_pdb in pdb_it0:
        iteration = 'it0'
        # path_pdb = 'run-true-interface/structures/it0/1GSL_1.pdb'
        cmd = f'python3 /trinity/login/apeliss/glycan-docking/scripts/capri-eval-fnat-corrected.py {ref_pdb} {path_pdb} --atoms-target C1,O4,C4,C5'
        profit_output = subprocess.getoutput(cmd)

        # I-RMSD
        i_rmsd_values = [line for line in profit_output.split(os.linesep)
                         if 'I-RMSD' in line]
        # we need the I-RMSD value
        try:
            i_rmsd = float(i_rmsd_values[-1].split()[-1])
        except IndexError:
            i_rmsd = 'NaN'

        # FNAT
        fnat_values = [line for line in profit_output.split(os.linesep)
                       if 'FNAT' in line]
        # we need the FNAT value
        try:
            fnat = float(fnat_values[-1].split()[-1])
        except IndexError:
            fnat = 'NaN'

        # L-RMSD
        l_rmsd_values = [line for line in profit_output.split(os.linesep)
                         if ' L-RMSD' in line]
        # we need the L-RMSD value
        try:
            l_rmsd = float(l_rmsd_values[-1].split()[-1])
        except IndexError:
            l_rmsd = 'NaN'

        # I-L-RMSD
        i_l_rmsd_values = [line for line in profit_output.split(os.linesep)
                           if 'I-L-RMSD' in line]
        # we need the I-L-RMSD value
        try:
            i_l_rmsd = float(i_l_rmsd_values[-1].split()[-1])
        except IndexError:
            i_l_rmsd = 'NaN'

        # ==========
        pdb = path_pdb.split('/')[-1]
        # pdb = '1GSL_1.pdb'
        rank = int(score_ranking_it0[pdb]['ranking'])
        haddock_score = float(score_ranking_it0[pdb]['haddock-score'])
        # ==========
        row = pdb_id, iteration, rank, i_rmsd, fnat, l_rmsd, i_l_rmsd, haddock_score, path_pdb
        rmsd_it0_data.append(row)

    with open('rmsd.csv', 'w', newline='') as csvf:
        rmsd_file = csv.writer(csvf)
        for line in rmsd_it0_data:
            rmsd_file.writerow(line)

    rmsd_it1_data = []
    for path_pdb in pdb_it1:
        iteration = 'it1'
        # path_pdb = 'run-true-interface/structures/it1/1GSL_1.pdb'
        cmd = f'python3 /trinity/login/apeliss/glycan-docking/scripts/capri-eval-fnat-corrected.py {ref_pdb} {path_pdb} --atoms-target C1,O4,C4,C5'
        profit_output = subprocess.getoutput(cmd)

        # I-RMSD
        i_rmsd_values = [line for line in profit_output.split(os.linesep)
                         if 'I-RMSD' in line]
        # we need the I-RMSD value
        try:
            i_rmsd = float(i_rmsd_values[-1].split()[-1])
        except IndexError:
            i_rmsd = 'NaN'

        # FNAT
        fnat_values = [line for line in profit_output.split(os.linesep)
                       if 'FNAT' in line]
        # we need the FNAT value
        try:
            fnat = float(fnat_values[-1].split()[-1])
        except IndexError:
            fnat = 'NaN'

        # L-RMSD
        l_rmsd_values = [line for line in profit_output.split(os.linesep)
                         if ' L-RMSD' in line]
        # we need the L-RMSD value
        try:
            l_rmsd = float(l_rmsd_values[-1].split()[-1])
        except:
            l_rmsd = 'NaN'

        # I-L-RMSD
        i_l_rmsd_values = [line for line in profit_output.split(os.linesep)
                           if 'I-L-RMSD' in line]
        # we need the I-L-RMSD value
        try:
            i_l_rmsd = float(i_l_rmsd_values[-1].split()[-1])
        except IndexError:
            i_l_rmsd = 'NaN'

        # ==========
        pdb = path_pdb.split('/')[-1]
        # pdb = '1GSL_1.pdb'
        rank = int(score_ranking_it1[pdb]['ranking'])
        haddock_score = float(score_ranking_it1[pdb]['haddock-score'])
        # ==========
        row = pdb_id, iteration, rank, i_rmsd, fnat, l_rmsd, i_l_rmsd, haddock_score, path_pdb
        rmsd_it1_data.append(row)

    with open('rmsd.csv', 'a', newline='') as csvf:
        rmsd_file = csv.writer(csvf)
        for line in rmsd_it1_data:
            rmsd_file.writerow(line)

    rmsd_itw_data = []
    for path_pdb in pdb_itw:
        iteration = 'itw'
        # path_pdb = 'run-true-interface/structures/itw/water/1GSL_1w.pdb'
        cmd = f'python3 /trinity/login/apeliss/glycan-docking/scripts/capri-eval-fnat-corrected.py {ref_pdb} {path_pdb} --atoms-target C1,O4,C4,C5'
        profit_output = subprocess.getoutput(cmd)

        # I-RMSD
        i_rmsd_values = [line for line in profit_output.split(os.linesep)
                         if 'I-RMSD' in line]
        # we need the I-RMSD value
        try:
            i_rmsd = float(i_rmsd_values[-1].split()[-1])
        except IndexError:
            i_rmsd = 'NaN'

        # FNAT
        fnat_values = [line for line in profit_output.split(os.linesep)
                       if 'FNAT' in line]
        # we need the FNAT value
        try:
            fnat = float(fnat_values[-1].split()[-1])
        except IndexError:
            fnat = 'NaN'

        # L-RMSD
        l_rmsd_values = [line for line in profit_output.split(os.linesep)
                         if ' L-RMSD' in line]
        # we need the L-RMSD value
        try:
            l_rmsd = float(l_rmsd_values[-1].split()[-1])
        except IndexError:
            l_rmsd = 'NaN'

        # I-L-RMSD
        i_l_rmsd_values = [line for line in profit_output.split(os.linesep)
                           if 'I-L-RMSD' in line]
        # we need the I-L-RMSD value
        try:
            i_l_rmsd = float(i_l_rmsd_values[-1].split()[-1])
        except IndexError:
            i_l_rmsd = 'NaN'

        # ==========
        pdb = path_pdb.split('/')[-1]
        # pdb = '1GSL_1w.pdb'
        rank = int(score_ranking_itw[pdb]['ranking'])
        haddock_score = float(score_ranking_itw[pdb]['haddock-score'])
        # ==========
        row = pdb_id, iteration, rank, i_rmsd, fnat, l_rmsd, i_l_rmsd, haddock_score, path_pdb
        rmsd_itw_data.append(row)

    with open('rmsd.csv', 'a', newline='') as csvf:
        rmsd_file = csv.writer(csvf)
        for line in rmsd_itw_data:
            rmsd_file.writerow(line)

# os.chdir(curr_dir)
# # fieldnames = ['pdb_id', 'iteration', 'rank', 'i-rmsd', 'fnat', 'l-rmsd', 'i-l-rmsd', 'haddock_score', 'path_to_the_structure']
# fieldnames = 'pdb_id,iteration,rank,i-rmsd,fnat,l-rmsd,i-l-rmsd,haddock_score,path_to_the_structure'
# with open('rmsd_all', 'w') as rmsd_all:
#     # rmsd_all = csv.writer(rmsd_all)
#     rmsd_all.write(fieldnames + '\n')
#     for pdb_id in directories_list:
#         w_direct = curr_dir + '/' + pdb_id + '/'
#         os.chdir(w_direct)
#         with open('rmsd.csv', 'r') as rmsd_ind:
#             for line in rmsd_ind:
#                 rmsd_all.write(line)
