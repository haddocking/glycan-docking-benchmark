# Calculate the pairwise distances between all atoms of the receptor and the ligand of a complex pdb file
import os
from math import sqrt

def interface_residues(r_dataset_file, l_dataset_file):  
    atdic_r = {}
    atdic_l = {}
    with open(r_dataset_file, "r") as r_fh:
        for line in r_fh.readlines():
            if line.startswith('ATOM'):
                resnum = int(line[23:26])
                atnum = int(line[7:11])
                x = float(line[31:38])
                y = float(line[39:46])
                z = float(line[47:54])
                atdic_r[atnum] = x, y, z, resnum
    with open(l_dataset_file, "r") as l_fh:
        for line in l_fh.readlines(): 
            if line.startswith('HETATM'):
                resnum = int(line[23:26])
                atnum = int(line[7:11])
                x = float(line[31:38])
                y = float(line[39:46])
                z = float(line[47:54])
                atdic_l[atnum] = x, y, z, resnum
    
    interface_r_list = []
    interface_l_list = []
    # Calculate the distance with the formula
    for at_r in atdic_r:
        xr, yr, zr = atdic_r[at_r][0], atdic_r[at_r][1], atdic_r[at_r][2]
        for at_l in atdic_l:
            xl, yl, zl = atdic_l[at_l][0], atdic_l[at_l][1], atdic_l[at_l][2]
            distance = sqrt((xr - xl)**2 + (yr - yl)**2 + (zr - zl)**2)
            if distance <= 3.9:
                interface_r_list.append(atdic_r[at_r][3])   # Append the resnum
                interface_l_list.append(atdic_l[at_l][3])
    
    # Remove duplicates and write to text files
    interface_r_list = list(set(interface_r_list))
    interface_l_list = list(set(interface_l_list))
    
    with open('receptor_con.txt', 'w') as r_con:
        interface_r = map(str, interface_r_list)
        interface_r = ' '.join(interface_r)
        r_con.write(interface_r)
    with open('ligand_con.txt', 'w') as l_con:
        interface_l = map(str, interface_l_list)
        interface_l = ' '.join(interface_l)
        l_con.write(interface_l)
