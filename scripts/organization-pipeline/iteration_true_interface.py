# Convert a list of residues from the receptor and ligand residues list files to a true-interface 
# restraints file (.tbl).

# Execution: python iteration-folders.py location/of/the/folders-for-each-structure-with-name-PDBID
# e.g., python .\iteration_true-interface.py .\HADDOCK-ready\ 

import os
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("residuefiles_folders", help="Residue files folders location")
args = parser.parse_args()
os.chdir(args.residuefiles_folders)

def generate_restraint(interface_r, interface_l, segid_r='A', segid_l='B'):
    # Parameters
    # interface_r: List of residue numbers of the receptor
    # interface_l: List of residue numbers of the ligand
    tbl_string = ''
    for resi_r in interface_r:
        tbl_string += 'assign (resi {} and segid {})'.format(resi_r, segid_r) + os.linesep
        tbl_string += '(' + os.linesep
        c = 0
        for resi_l in interface_l:
            tbl_string += '       (resi {} and segid {})'.format(resi_l, segid_l) + os.linesep
            c += 1
            if c != len(interface_l):
                tbl_string += '        or' + os.linesep
        tbl_string += ') 2.0 2.0 0.0' + os.linesep

    for resi_l in interface_l:
        tbl_string += 'assign (resi {} and segid {})'.format(resi_l, segid_l) + os.linesep
        tbl_string += '(' + os.linesep
        c = 0
        for resi_r in interface_r:
            tbl_string += '       (resi {} and segid {})'.format(resi_r, segid_r) + os.linesep
            c += 1
            if c != len(interface_r):
                tbl_string += '        or' + os.linesep
        tbl_string += ') 2.0 2.0 0.0' + os.linesep

    return tbl_string


def load_file(residue_file):
    """Read a space-separated file and return a list of integers."""
    output_list = []
    with open(residue_file, 'r') as fh:
        for line in fh.readlines():
            # all the residues are in one single line
            residue_str_list = line.split()  # remove spaces
            residue_int_list = map(int, residue_str_list)  # convert to integers
            # add it to the residue_list
            for res in residue_int_list:
                output_list.append(res)
    return output_list


# def main():
    # parser = argparse.ArgumentParser(description='Generate a .tbl file based on two files containing the interface residues')
    # parser.add_argument('receptor_f', help='File containing the receptor residues separated by space')
    # parser.add_argument('ligand_f', help='File containing the ligand residues separated by space')
    # parser.add_argument('output_f', help='Output file that will contain the restraints')
    # args = parser.parse_args()

curr_dir = os.getcwd()
directories_list = []
for subdir, dirs, files in os.walk(curr_dir):
    for direct in dirs:
        directories_list.append(direct)

for structure_direct in directories_list:
    w_direct = curr_dir + '/' + structure_direct + '/'
    os.chdir(w_direct)
    # print(os.getcwd())
 
    files_list = []
    
    for filename in os.listdir():
        files_list.append(filename)
    # print(files_list)
    
    for receptor_con in files_list:
        if receptor_con.startswith("receptor_con"):
            for ligand_con in files_list:
                if ligand_con.startswith("ligand_con"):
                    # do some check before we start running

                    # check if the receptor file exist
                    if not os.path.isfile(receptor_con):
                        print(f'Receptor file {receptor_con} not found')
                        # parser.print_help()  # this will print the usage
                        # sys.exit()
                    
                    # check if the ligand file exist
                    if not os.path.isfile(ligand_con):
                        print(f'Ligand file {ligand_con} not found')
                    #     parser.print_help()
                    #     sys.exit()
                    
                    # check if this output file already exists
                    # if os.path.isfile(args.output_f):
                    #     print(f'Output file {args.output_f} already exists!')
                    #     sys.exit()

                    interface_r = load_file(receptor_con)
                    interface_l = load_file(ligand_con)

                    restraint_str = generate_restraint(interface_r, interface_l)

                    with open("true_interface.tbl", 'w') as out_fh:
                        out_fh.write(restraint_str)
                    out_fh.close()


# if __name__ == '__main__':
#     main()
