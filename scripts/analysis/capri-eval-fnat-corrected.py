from _operator import itemgetter
import argparse
import shutil
import subprocess
import glob
import os
import sys
import random
import logging
from itertools import groupby
# import operator


def split_chain(pdbf):
    # TODO: Properly implement the pdb-tools here
    """Split PDB by chain."""
    logging.debug(f'Splitting {pdbf} by ChainID')
    prev_chain, chain_ids, chain_atoms = None, [], {}
    cur_chain = None
    for line in open(pdbf):
        if 'ATOM' in line[:4]:
            if prev_chain != line[21]:
                if not line[21] in chain_atoms:
                    cur_chain = chain_atoms[line[21]] = []
                else:
                    cur_chain = chain_atoms[line[21]]
                cur_chain.append(line)
                prev_chain = line[21]
                chain_ids.append(line[21])
            else:
                cur_chain.append(line)

    # Output chains to files
    pdb_dic = {}
    for c_id in chain_ids:
        name = pdbf.split('.pdb')[0] + '_' + c_id + '.pdb'
        pdb_dic[c_id] = name
        out = open(name, 'w')
        out.write(''.join(chain_atoms[c_id]))
        out.write('END\n')
        out.close()

    return pdb_dic


def load_seq(prot_fname):
    """Read a PDB file and return the sequence."""
    logging.debug(f'Retrieving sequence of {prot_fname}')
    aa_dic = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
              'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
              'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
              'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
              'DC': 'C', 'DA': 'A', 'DG': 'G', 'DT': 'T', 'ADE': 'A',
              'THY': 'T', 'GUA': 'G', 'CYT': 'C'}
    seq_dic = {}
    with open(prot_fname, 'r') as fh:
        for line in fh.readlines():
            if 'ATOM' in line[:4]:
                chain = line[21]
                resnum = int(line[22:26])
                resname = line[17:20].split()[0]
                try:
                    _ = seq_dic[chain]
                except KeyError:
                    seq_dic[chain] = {}
                try:
                    name = aa_dic[resname]
                except KeyError:
                    name = 'X'
                seq_dic[chain][resnum] = name
    return seq_dic


def align(prota_fname, protb_fname):
    """Structurally align proteins and return a matching numbering dictionary.

    {'A': {protA_res_1: protB_res_42, protA_res_2: protB_res_43, etc }}
    """
    logging.debug('Aligning structures with LOVOALIGN')
    # FIXME: Why splitting here?
    protein_a_dic = split_chain(prota_fname)
    protein_b_dic = split_chain(protb_fname)

    # check if chain ids match
    if protein_a_dic.keys() != protein_b_dic.keys():
        _msg = f'ChainIDs of {prota_fname} and {protb_fname} do not match!'
        logging.error(_msg)
        sys.exit()

    numbering_dic = {}
    for chain in protein_a_dic.keys():
        logging.debug(f'Structurally aligning chain {chain}')
        numbering_dic[chain] = {}
        cmd = (f'lovoalign -p1 {protein_a_dic[chain]} '
               f'-p2 {protein_b_dic[chain]} '
               f'-c1 {chain} -c2 {chain}')

        logging.debug(f'Command is: {cmd}')
        lovoalign_output = subprocess.getoutput(cmd).split(os.linesep)

        # find out where the alignment starts and ends
        alignment_pass = True
        for i, line in enumerate(lovoalign_output):
            if 'SEQUENCE ALIGNMENT' in line:
                # there are 2 extra white lines after this header
                alignment_start_index = i + 2
            elif 'FINAL' in line:
                # there are 2 extra white lines after this header
                alignment_end_index = i - 2
            elif 'ERROR' in line:
                failed_pdb = line.split()[-1]
                _msg = (f'LovoAlign could not read {failed_pdb} '
                        'is it a ligand?')
                logging.warning(_msg)
                alignment_pass = False
                pa_seqdic = load_seq(prota_fname)
                for elem in [k for k in pa_seqdic[chain]]:
                    numbering_dic[chain][elem] = elem

        if not alignment_pass:
            # This alignment failed, move on to the next
            logging.warning(f'Skipping alignment of chain {chain}')
            continue

        aln_l = lovoalign_output[alignment_start_index:alignment_end_index]

        # remove the line between the alignment segments
        alignment = [aln_l[i:i + 3][:2] for i in range(0, len(aln_l), 3)]

        logging.debug('Reading alignment and matching numbering')
        for element in alignment:
            line_a, line_b = element

            resnum_a, seq_a, _ = line_a.split()
            resnum_b, seq_b, _ = line_b.split()

            resnum_a = int(resnum_a) - 1
            resnum_b = int(resnum_b) - 1

            for resname_a, resname_b in zip(seq_a, seq_b):
                if resname_a != '-':
                    resnum_a += 1

                if resname_b != '-':
                    resnum_b += 1

                if resname_a != '-' and resname_b != '-':
                    numbering_dic[chain][resnum_b] = resnum_a

    return numbering_dic


def output_renumbered(prot, numbering_dic):
    #
    renumb_pdb_l = []
    for line in open(prot):
        if 'ATOM' in line[:4]:
            chain = line[21]
            resnum = int(line[22:26])
            try:
                new_res = numbering_dic[chain][resnum]
                n_l = line[:22] + f'{new_res:>4}' + line[26:]
                renumb_pdb_l.append(n_l)
            except KeyError:
                print(f'Residue {resnum} does not have a match'
                      ', it will be skipped')
                pass


    outputf = prot.replace('.pdb', '-renum.pdb')
    out = open(outputf, 'w')
    out.write(''.join(renumb_pdb_l))
    out.close()
    print(prot)


def run_contacts(pdbf, cutoff):
    """Run the contacts script."""
    cmd = f'contact {pdbf} {cutoff}'
    out = subprocess.getoutput(cmd).split(os.linesep)
    return out


def identify_inteface(pdbf, cutoff):
    """Identify the interfacing residues given a cutoff."""
    contacts_l = run_contacts(pdbf, cutoff)

    interface_dic = {}
    for line in contacts_l:
        data = line.split()
        resnum_a = int(data[0])
        chain_a = data[1]
        # atom_a = data[2]
        resnum_b = int(data[3])
        chain_b = data[4]
        # atom_b = data[5]
        distance = float(data[6])

        # One way
        try:
            _ = interface_dic[chain_a]
        except KeyError:
            interface_dic[chain_a] = {}
        try:
            _ = interface_dic[chain_a][chain_b]
        except KeyError:
            interface_dic[chain_a][chain_b] = []

        if distance <= cutoff:
            if resnum_a not in interface_dic[chain_a][chain_b]:
                interface_dic[chain_a][chain_b].append(resnum_a)

        # other way
        try:
            _ = interface_dic[chain_b]
        except KeyError:
            interface_dic[chain_b] = {}
        try:
            _ = interface_dic[chain_b][chain_a]
        except KeyError:
            interface_dic[chain_b][chain_a] = []

        if float(distance) <= cutoff:
            if resnum_b not in interface_dic[chain_b][chain_a]:
                interface_dic[chain_b][chain_a].append(resnum_b)

    # sort residue lists
    ninterface_dic = {}
    for a in interface_dic:
        ninterface_dic[a] = {}
        for b in interface_dic[a]:
            reslist = interface_dic[a][b]
            reslist.sort()
            ninterface_dic[a][b] = reslist

    return ninterface_dic


def get_range(data):
    ranges = []
    for k, g in groupby(enumerate(data), lambda x: x[0] - x[1]):
        group = (map(itemgetter(1), g))
        group = list(map(int, group))
        ranges.append((group[0], group[-1]))
    return ranges


def retrieve_izone(c_dic, numbering_dic):
    # based on the reference interface, create izone
    izone_l = []
    for chain in c_dic:
        ref_dic = {}
        for bound_res in list(c_dic[chain].items())[0][1]:
            try:
                ub = numbering_dic[chain][bound_res]
                ref_dic[bound_res] = ub
            except KeyError:
                pass

        for bound_range in get_range(ref_dic.keys()):
            unbound_res_l = []
            for bound_res in range(bound_range[0], bound_range[1] + 1):
                unbound_res_l.append(ref_dic[bound_res])

            for unbound_range in get_range(unbound_res_l):
                bound_res_l = []
                for unbound_res in range(unbound_range[0], unbound_range[1] + 1):
                    bound_res_l.append(list(ref_dic.keys())[list(ref_dic.values()).index(unbound_res)])

                range_a = get_range(bound_res_l)[0]  # bound
                range_b = unbound_range

                izone_str = ('ZONE '
                             f'{chain}{range_a[0]}-{chain}{range_a[1]}:'
                             f'{chain}{range_b[0]}-{chain}{range_b[1]}')
                izone_l.append(izone_str)

    return izone_l


def run_profit(cmd):
    return subprocess.getoutput(f'echo "{cmd}" | profit').split(os.linesep)


def calc_i_rmsd(receptor_mol, ligand_mol, receptor_atoms, ligand_atoms,
                numbering_dic, cutoff=10.0):

    contact_dic_a = identify_inteface(receptor_mol, cutoff=cutoff)
    izone_l = retrieve_izone(contact_dic_a, numbering_dic)
    izone_str = os.linesep.join(izone_l)

    cmd = f'REFE {receptor_mol}' + os.linesep
    cmd += f'MOBI {ligand_mol}' + os.linesep
    cmd += f'ATOMS {receptor_atoms},{ligand_atoms}' + os.linesep
    cmd += 'ZONE CLEAR' + os.linesep
    cmd += f'{izone_str}' + os.linesep
    # cmd += 'STATUS' + os.linesep
    cmd += 'FIT' + os.linesep
    cmd += 'QUIT' + os.linesep

    with open('i-rmsd.dbg', 'w') as fh:
        fh.write(cmd)

    with open('izone', 'w') as fh:
        fh.write(izone_str)

    profit_output = run_profit(cmd)
    with open('i-rmsd.out', 'w') as fh:
        for line in profit_output:
            fh.write(line + os.linesep)
            if 'Error' in line:
                _msg = 'PROFIT raised an error! Check i-rmsd.out'
                logging.warning(_msg)

    irmsd = None
    try:
        irmsd = float(profit_output[-2].split()[-1])
    except KeyError:
        _msg = 'Something went wrong when running PROFIT, check i-rmsd.dbg'
        logging.error(_msg)
        sys.exit()

    return irmsd


def calc_fnat(receptor_mol, ligand_mol, numbering_dic, cutoff=5.0):
    """Calculate the frequency of native contacts."""
    receptor_contacts = run_contacts(receptor_mol, cutoff=cutoff)
    ligand_contacts = run_contacts(ligand_mol, cutoff=cutoff)

    receptor_con_list = []
    ligand_con_list = []
    for line in receptor_contacts:
        resnum_x, chain_x, _, resnum_y, chain_y, _, _ = line.split()
        try:
            resnum_x = str(numbering_dic[chain_x][int(resnum_x)])
            resnum_y = str(numbering_dic[chain_y][int(resnum_y)])
        except KeyError:
            # one of the residues present in this contact was not matched to
            #  the target
            continue

        receptor_con_list.append((resnum_x, resnum_y))
    receptor_con_list = set(receptor_con_list)

    for line in ligand_contacts:
        try:
            resnum_x, _, _, resnum_y, _, _, _ = line.split()

            ligand_con_list.append((resnum_x, resnum_y))
        except ValueError:
            pass
    ligand_con_list = set(ligand_con_list)

    try:
        fnat = (len(receptor_con_list & ligand_con_list) /
                len(receptor_con_list))
    except ZeroDivisionError:
        # No contacts were matched
        fnat = .0

    return fnat


def calc_l_rmsd(receptor_mol, ligand_mol, receptor_atoms, ligand_atoms,
                numbering_dic):
    """Calculate the ligand-RMSD."""
    # This is done by aligning on the receptor and calculating over the atoms
    #  of the ligand
    lrmsd = None

    chain_l = list(numbering_dic.keys())
    chain_l.sort()
    # Warning, receptor will always be the first chain!
    receptor_chain = chain_l[0]
    zone_dic = {}
    # Create the lzone
    for chain in numbering_dic:
        bound_reslist = list(numbering_dic[chain].keys())
        unbound_reslist = list(numbering_dic[chain].values())

        zone_dic[chain] = []
        # for each bound residue range
        for bound_range in get_range(numbering_dic[chain]):

            # found the unbound equivalents
            unbound_res_l = []
            for bound_res in range(bound_range[0], bound_range[1] + 1):
                unbound_res = numbering_dic[chain][bound_res]
                unbound_res_l.append(unbound_res)

            # do the other way around (?)
            #  based in the unbound range, find the bound equivalent
            #   (I could not figure a way to do this in any other way,
            #    but for sure it could be improved)
            for unbound_rng in get_range(unbound_res_l):
                bound_res_l = []
                for unbound_res in range(unbound_rng[0], unbound_rng[1] + 1):
                    unbound_index = unbound_reslist.index(unbound_res)
                    bound_res_l.append(bound_reslist[unbound_index])

                range_a = get_range(bound_res_l)[0]  # bound
                range_b = unbound_rng

                receptor_zone_str = (f'ZONE '
                                     f'{chain}{range_a[0]}-{chain}{range_a[1]}'
                                     ':'
                                     f'{chain}{range_b[0]}-{chain}{range_b[1]}'
                                     )
                zone_dic[chain].append(receptor_zone_str)

    receptor_zone = os.linesep.join(zone_dic[receptor_chain])
    lzone = ''
    cmd = f'REFE {receptor_mol}' + os.linesep
    cmd += f'MOBI {ligand_mol}' + os.linesep
    cmd += f'ATOMS {receptor_atoms}' + os.linesep
    cmd += receptor_zone + os.linesep
    cmd += 'FIT' + os.linesep

    # cmd += 'RATOMS ^H*' + os.linesep
    cmd += f"RATOMS {ligand_atoms}" + os.linesep

    for ligand_chain in zone_dic:
        if ligand_chain != receptor_chain:
            for ligand_zone in zone_dic[ligand_chain]:
                cmd += f'R{ligand_zone}' + os.linesep
                lzone += f'R{ligand_zone}' + os.linesep
    cmd += 'ZONE CLEAR' + os.linesep
    cmd += 'QUIT'

    with open('lrmsd.dbg', 'w') as fh:
        fh.write(cmd)

    with open('lzone', 'w') as fh:
        fh.write(lzone)

    profit_output = run_profit(cmd)
    with open('l-rmsd.out', 'w') as fh:
        for line in profit_output:
            fh.write(line + os.linesep)
            if 'Error' in line:
                _msg = 'PROFIT raised an error! Check l-rmsd.out'
                logging.warning(_msg)
    try:
        lrmsd = float([e for e in profit_output if 'RMS' in e][-1].split()[-1])
    except KeyError:
        _msg = 'Something went wrong when running PROFIT, check l-rmsd.out'
        logging.error(_msg)
        sys.exit()

    return lrmsd


def calc_i_l_rmsd(receptor_mol, ligand_mol, receptor_atoms, ligand_atoms,
                  numbering_dic, cutoff=5.0):
    """Calculate the interface-ligand-RMSD."""
    # This is done by aligning on the receptor interface and calculating
    # over the atoms of the ligand
    chain_l = list(numbering_dic.keys())
    chain_l.sort()
    # Warning, receptor will always be the first chain!
    receptor_chain = chain_l[0]
    ligand_chain = chain_l[1]
    zone_dic = {}
    # Create the lzone
    for chain in numbering_dic:
        bound_reslist = list(numbering_dic[chain].keys())
        unbound_reslist = list(numbering_dic[chain].values())

        zone_dic[chain] = []
        # for each bound residue range
        for bound_range in get_range(numbering_dic[chain]):

            # found the unbound equivalents
            unbound_res_l = []
            for bound_res in range(bound_range[0], bound_range[1] + 1):
                unbound_res = numbering_dic[chain][bound_res]
                unbound_res_l.append(unbound_res)

            # do the other way around (?)
            #  based in the unbound range, find the bound equivalent
            #   (I could not figure a way to do this in any other way,
            #    but for sure it could be improved)
            for unbound_rng in get_range(unbound_res_l):
                bound_res_l = []
                for unbound_res in range(unbound_rng[0], unbound_rng[1] + 1):
                    unbound_index = unbound_reslist.index(unbound_res)
                    bound_res_l.append(bound_reslist[unbound_index])

                range_a = get_range(bound_res_l)[0]  # bound
                range_b = unbound_rng

                receptor_zone_str = (f'ZONE '
                                     f'{chain}{range_a[0]}-{chain}{range_a[1]}'
                                     ':'
                                     f'{chain}{range_b[0]}-{chain}{range_b[1]}'
                                     )
                zone_dic[chain].append(receptor_zone_str)

    contact_dic_a = identify_inteface(receptor_mol, cutoff=cutoff)
    izone_l = retrieve_izone(contact_dic_a, numbering_dic)
    izone_str = os.linesep.join([j for j in izone_l if ligand_chain not in j])

    cmd = f'REFE {receptor_mol}' + os.linesep
    cmd += f'MOBI {ligand_mol}' + os.linesep
    cmd += f'ATOMS {receptor_atoms}' + os.linesep
    cmd += 'ZONE CLEAR' + os.linesep
    cmd += f'{izone_str}' + os.linesep
    # cmd += 'STATUS' + os.linesep
    cmd += 'FIT' + os.linesep

    cmd += f"RATOMS {ligand_atoms}" + os.linesep

    lzone = ''
    for ligand_chain in zone_dic:
        if ligand_chain != receptor_chain:
            for ligand_zone in zone_dic[ligand_chain]:
                cmd += f'R{ligand_zone}' + os.linesep
                lzone += f'R{ligand_zone}' + os.linesep
    cmd += 'QUIT' + os.linesep

    with open('i-l-rmsd.dbg', 'w') as fh:
        fh.write(cmd)

    with open('lzone', 'w') as fh:
        fh.write(lzone)

    profit_output = run_profit(cmd)
    with open('i-l-rmsd.out', 'w') as fh:
        for line in profit_output:
            fh.write(line + os.linesep)
            if 'Error' in line:
                _msg = 'PROFIT raised an error! Check i-l-rmsd.out'
                logging.warning(_msg)
    try:
        lirmsd = float([e for e in profit_output if 'RMS' in e][-1].split()[-1])
    except KeyError:
        _msg = 'Something went wrong when running PROFIT, check i-l-lrmsd.out'
        logging.error(_msg)
        sys.exit()

    return lirmsd


def clean(prot_a, prot_b):
    """Remove files generated by this script."""
    pdb_l = [glob.glob(f"{e.split('.pdb')[0]}_*") for e in [prot_a, prot_b]]
    pdb_l = [x for xs in pdb_l for x in xs] + glob.glob('*flatnum*')
    for p in set(pdb_l):
        os.system(f'rm {p}')


def scramble_prot(prot):
    """Scramble protein numbering for debug purposes."""
    seq = load_seq(prot)
    chain_resdic = dict([(c, list(seq[c].keys())) for c in seq])
    new_resdic = {}
    for chain in chain_resdic:
        new_resdic[chain] = {}
        reslist = chain_resdic[chain]
        shuffled_reslist = reslist.copy()
        random.shuffle(shuffled_reslist)
        new_resdic[chain] = dict(zip(reslist, shuffled_reslist))

    scrambled_prot = []
    for line in open(prot):
        if 'ATOM' in line[:4]:
            chain = line[21]
            resnum = int(line[22:26])
            new_resnum = new_resdic[chain][resnum]
            n_l = line[:22] + '{0:>4}'.format(new_resnum) + line[26:]
            scrambled_prot.append(n_l)

    outname = '{}_scramb.pdb'.format(prot.split('.pdb')[0])
    out = open(outname, 'w')
    out.write(''.join(scrambled_prot))
    out.close()

    return outname


# TODO: Use pdb-tools instead
# def flatten_numbers(prot):
#     """Flatten number insertions."""
#     # look for residue insertions
#     # 10, <10A, 10B, 10C>, 12
#     # create a numbering dictionary taking into account insertions
#     resdic = {}
#     chain_l = []
#     incr = None
#     for l in open(prot):
#         if 'ATOM' in l[:4]:
#             chain = l[21]
#             if not chain in chain_l:
#                 incr = 0
#                 chain_l.append(chain)
#             try:
#                 _ = resdic[chain]
#             except KeyError:
#                 resdic[chain] = {}
#             ori_resnum = int(l[22:26])
#             icode = l[26]
#             if not icode.isspace():
#                 incr += 1
#             resnum = ori_resnum + incr
#             resdic[chain][(ori_resnum, icode)] = resnum
#     flatf_l = []
#     for l in open(prot):
#         if 'ATOM' in l[:4]:
#             chain = l[21]
#             resnum = int(l[22:26])
#             icode = l[26]
#             new_res = resdic[chain][(resnum, icode)]
#             n_l = l[:22] + '{:>4}'.format(new_res) + ' ' + l[27:]
#             flatf_l.append(n_l)
#     outputf = prot.split('.pdb')[0] + '_flatnum.pdb'
#     out = open(outputf,'w')
#     out.write(''.join(flatf_l))
#     out.close()
#     return outputf


def calc_dockq(fnat, irms, lrms):
    """Calculate the DockQ metric."""
    # This equation is based on:
    # https://github.com/bjornwallner/DockQ/blob/3735c160050f1e9128d2ccb23a0a1945aa98b5b2/DockQ.py#L361
    dockq = (fnat + 1 / (1 + (irms / 1.5) * (irms / 1.5)) +
             1 / (1 + (lrms / 8.5) * (lrms / 8.5))) / 3
    return dockq


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("reference", help="Reference complex")
    parser.add_argument("target", help="Target complex")

    parser.add_argument("--atoms-reference",
                        help="Atom selection for the reference complex",
                        type=str,
                        default='C,N,CA,O')
    parser.add_argument("--atoms-target",
                        help="Atom selection for the target complex",
                        type=str,
                        default='C,N,CA,O')

    parser.add_argument("--flat",
                        help=("Flatten numbering, use this if your complex has"
                              " insertions ex: 10, 10A, 10B, 11"))

    levels = ('DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL')
    parser.add_argument('--log-level', default='INFO', choices=levels)
    args = parser.parse_args()

    logging.basicConfig(level=args.log_level,
                        format=('[%(asctime)s] %(funcName)s():L%(lineno)d '
                                '%(levelname)s - %(message)s'),
                        datefmt='%d/%m/%Y %H:%M:%S')

    error_check = False
    for exec in ['profit', 'contact', 'lovoalign']:
        if not shutil.which(exec):
            logging.error(f' {exec} not found in $PATH')
            error_check = True
    if error_check:
        sys.exit()

    protein_a_fname = args.reference
    protein_b_fname = args.target

    atoms_a = args.atoms_reference
    atoms_b = args.atoms_target

    logging.info(f'Receptor: {protein_a_fname} Atoms: {atoms_a}')
    logging.info(f'Ligand: {protein_b_fname} Atoms: {atoms_b}')

    # if args.flat:
    #     protein_a_fname = flatten_numbers(protein_a_fname)
    #     protein_b_fname = flatten_numbers(protein_b_fname)

    # DEBUG ONLY!
    # pb = scramble_prot(pb)

    num_dic = align(protein_a_fname, protein_b_fname)

    # renumber the reference based on the target
    #output_renumbered(protein_b_fname, num_dic)
    # sys.exit()

    i_rmsd = calc_i_rmsd(protein_a_fname,
                         protein_b_fname,
                         atoms_a,
                         atoms_b,
                         num_dic,
                         cutoff=6.0)

    fnat = calc_fnat(protein_a_fname,
                     protein_b_fname,
                     num_dic,
                     cutoff=4.0)

    l_rmsd = calc_l_rmsd(protein_a_fname,
                         protein_b_fname,
                         atoms_a,
                         atoms_b,
                         num_dic)

    i_l_rmsd = calc_i_l_rmsd(protein_a_fname,
                             protein_b_fname,
                             atoms_a,
                             atoms_b,
                             num_dic,
                             cutoff=6.0)

    dockq = calc_dockq(i_rmsd, fnat, l_rmsd)

    logging.info(f'I-RMSD: {i_rmsd:.2f}')
    logging.info(f'FNAT: {fnat:.2f}')
    logging.info(f'L-RMSD: {l_rmsd:.2f}')
    logging.info(f'I-L-RMSD: {i_l_rmsd:.2f}')
    logging.info(f'DockQ: {dockq:.2f} (untested)')

    clean(protein_a_fname, protein_b_fname)
