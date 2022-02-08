from selenium import webdriver
import pandas as pd
import pathlib
import tarfile
import time
import os
import re


if __name__ == "__main__":

    chrome_options = webdriver.ChromeOptions()
    prefs = {"profile.managed_default_content_settings.images": 2,
             "download.default_directory": str(pathlib.Path.cwd())}
    chrome_options.add_experimental_option("prefs", prefs)

    driver = webdriver.Chrome('C:\\Users\\Angela\\chromedriver_win32\\'
                              'chromedriver.exe')
    dataset = pd.read_excel('C:\\Users\\Angela\\OneDrive - Universiteit '
                            'Utrecht\\Major Research Project\\glycan-benchmark'
                            '\\01-organization\\unbound\\gd-dataset.xlsx')
    dataset.to_csv('dataset.csv', index=None, header=True)
    pd.set_option('max_colwidth', 1000)

    glycans = {}
    for pdbid in dataset['PDB_ID']:
        sequence = str(dataset[dataset['PDB_ID'] == pdbid]['Glycoligand_Sequence'])
        sequence = sequence[5:].strip()
        sequence = sequence.partition('\n')
        glycans[pdbid] = sequence[0]

    for pdbid, sequence in glycans.items():
        glycam_str = sequence.replace('>', '')
        url = "http://glycam.org/url?glycam=" + glycam_str + '-OH'
        driver.get(url)
        driver.get('http://glycam.org/tools/molecular-dynamics/oligosaccharide'
                   '-builder/download-files')

        wait = True
        while wait:
            time.sleep(15)
            try:
                driver.find_element_by_xpath("//input[@type='submit' and"
                                             " @value='Download All Stru"
                                             "ctures']").click()
                wait = False
            except:
                pass

            # open file
            time.sleep(20)
            file = tarfile.open('C:\\Users\\Angela\\Downloads\\'
                                'structures.tar.gz')
            # extracting file
            file.extractall('C:\\Users\\Angela\\OneDrive - Universiteit '
                            'Utrecht\\Major Research Project\\glycan-benchmark'
                            '\\01-organization\\unbound\\glycan-strct\\' +
                            pdbid)
            downloaded_files = os.listdir('C:\\Users\\Angela\\OneDrive - '
                                          'Universiteit Utrecht\\Major '
                                          'Research Project\\glycan-benchmark'
                                          '\\01-organization\\unbound\\glycan-'
                                          'strct\\' + pdbid)
            structures = []
            for pdb_file in downloaded_files:
                file_name = re.compile(r'\d.pdb')
                matches = file_name.findall(pdb_file)
                if matches:
                    structures.append(pdb_file)

            for structure in structures:
                if len(structures) == 1:
                    with open('C:\\Users\\Angela\\OneDrive - Universiteit'
                              ' Utrecht\\Major Research Project\\glycan-'
                              'benchmark\\01-organization\\unbound\\glycan'
                              '-strct\\' + pdbid + '\\' + pdbid +
                              '_l_u.pdb', 'w') as clean_ligand:
                        with open('C:\\Users\\Angela\\OneDrive - '
                                  'Universiteit Utrecht\\Major Research '
                                  'Project\\glycan-benchmark\\01-'
                                  'organization\\unbound\\glycan-strct\\' +
                                  pdbid + '\\' + structure, 'rt') as downloaded_ligand:
                            conversion = {
                                    'GB': 'BGC',
                                    'MB': 'BMA',
                                    'FA': 'FCA',
                                    'fA': 'FUC',
                                    'FB': 'FCB',
                                    'LB': 'GAL',
                                    'LA': 'GLA',
                                    'GA': 'GLC',
                                    'MA': 'MAN',
                                    'YB': 'NAG',
                                    'YA': 'NDG',
                                    'VB': 'NGA',
                                    'SA': 'SIA',
                                    'SB': 'SIB',
                                    'XB': 'XYP'
                                    }
                            pattern = re.compile('^HETATM')
                            for line in downloaded_ligand:
                                if re.search(pattern, line):
                                    resid = str(line[17:20])
                                    if resid != 'ROH':
                                        for glycam_notation in conversion:
                                            if resid[1:] == glycam_notation:
                                                clean_ligand.write(line[0:17] + conversion[glycam_notation] + line[20:])
                                elif line.startswith('END'):
                                    clean_ligand.write(line)
                else:
                    with open('C:\\Users\\Angela\\OneDrive - Universiteit'
                              ' Utrecht\\Major Research Project\\glycan-'
                              'benchmark\\01-organization\\unbound\\glycan'
                              '-strct\\' + pdbid + '\\' + pdbid +
                              '_l_u_' + structure[0] + '.pdb', 'w') as clean_ligand:
                        with open('C:\\Users\\Angela\\OneDrive - '
                                  'Universiteit Utrecht\\Major Research '
                                  'Project\\glycan-benchmark\\01-'
                                  'organization\\unbound\\glycan-strct\\' +
                                  pdbid + '\\' + structure, 'rt') as downloaded_ligand:
                            conversion = {
                                    'GB': 'BGC',
                                    'MB': 'BMA',
                                    'FA': 'FCA',
                                    'fA': 'FUC',
                                    'FB': 'FCB',
                                    'LB': 'GAL',
                                    'LA': 'GLA',
                                    'GA': 'GLC',
                                    'MA': 'MAN',
                                    'YB': 'NAG',
                                    'YA': 'NDG',
                                    'VB': 'NGA',
                                    'SA': 'SIA',
                                    'SB': 'SIB',
                                    'XB': 'XYP'
                                    }
                            pattern = re.compile('^HETATM')
                            for line in downloaded_ligand:
                                if re.search(pattern, line):
                                    resid = str(line[17:20])
                                    if resid != 'ROH':
                                        for glycam_notation in conversion:
                                            if resid[1:] == glycam_notation:
                                                clean_ligand.write(line[0:17] + conversion[glycam_notation] + line[20:])
                                elif line.startswith('END'):
                                    clean_ligand.write(line)

            file.close()

            os.unlink('C:\\Users\\Angela\\Downloads\\structures.tar.gz')
