import Bio
from Bio.PDB import PDBParser
from Bio.PDB import PDBList
import urllib
import pandas as pd
import xml.etree.ElementTree as et
from urllib.request import urlopen
import requests
import os
import pymol
from pymol import cmd
import operator
from biopandas.pdb import PandasPdb


def download_pdb_chains(dali_search_output, base_dir):
    print("Start downloading - wait for finish signal\n\n")
    cTerminalDF = pd.read_csv(dali_search_output)
    cTerminalList = cTerminalDF["Chains"].tolist()
    space = " "
    command_pdb_sel_chain = "pdb_selchain"
    command_pdb_fetch = "pdb_fetch"
    redirect_output = ">"
    pdb_file_extension = ".pdb"
    slash = "/"
    relative_dir = "." + slash
    output_dir = relative_dir + base_dir + slash
    pdb_chain_downloaded_list = []
    for i in cTerminalList:
        chain = i[4:6]
        pdb_file = i[0:4]
        pdb_chain_downloaded = output_dir + pdb_file + chain + pdb_file_extension
        pdb_file_with_extension = output_dir + pdb_file + pdb_file_extension
        pdb_chain_downloaded_list.append(pdb_chain_downloaded)
        fetch_command = command_pdb_fetch + space + pdb_file + \
            space + redirect_output + space + pdb_file_with_extension
        chain_command = command_pdb_sel_chain + space + chain + space + \
            pdb_file_with_extension + space + redirect_output + space + pdb_chain_downloaded
        # FETCH THE PDB FILE
        os.system(fetch_command)
        # CREATE THE CHAIN ONLY PDB FILE
        os.system(chain_command)
    print("Done downloading\n\n")
    return pdb_chain_downloaded_list, cTerminalList, output_dir


cTerminal_chain_file_list, cTerminal_chain_name_list, cTerminal_output_dir = download_pdb_chains(
    "./CTerminalData.csv", "cTerminal-RNF41-Dali-Output-active-sites-pdb")
zf1_chain_file_list, zf1_chain_name_list, zf1_output_dir = download_pdb_chains(
    "./zf1.csv", "zf1-RNF41-Dali-Output-active-sites-pdb")
zf2_chain_file_list, zf2_chain_name_list, zf2_output_dir = download_pdb_chains(
    "./zf2.csv", "zf2-RNF41-Dali-Output-active-sites-pdb")


def align_calculate_RMSD(base_file_name, chain_file_list, output_dir, chain_name_list):
    pdb_file_extension = '.pdb'
    c_terminal_rnf41_pdb = output_dir + base_file_name + pdb_file_extension
    cmd.load(c_terminal_rnf41_pdb)
    chain_name_index = 0
    align_output_list = []
    for chain_file in chain_file_list:
        chain_name = chain_name_list[chain_name_index]
        chain_name_index += 1
        cmd.load(chain_file)
        output = cmd.align(base_file_name, chain_name, cycles=0, transform=0)
        align_output_list.append(output)
        cmd.delete(chain_file)
    cmd.delete(c_terminal_rnf41_pdb)
    return align_output_list


def print_align_output(output, chain_name, detailed):
    if (detailed):
        print("Align output for ", chain_name, ":")
        print("RMSD after refinement:", output[0])
        print("Number of aligned atoms after refinement:", output[1])
        print("Number of refinement cycles:", output[2])
        print("RMSD before refinement:", output[3])
        print("Number of aligned atoms before refinement:", output[4])
        print("Raw alignment score:", output[5])
        print("Number of residues aligned:", output[6])
    else:
        print("Align output for ", chain_name, ":", output)


def build_dictionary(output, chain_name, dict_name):
    dict_name.update({chain_name: output[0]})


def get_RSMD_sorted(dict_name):
    sorted_d = dict(sorted(dict_name.items(), key=operator.itemgetter(1)))
    return sorted_d


def get_alignment_RMSDs(input_protein, chain_file_list, output_dir, chain_name_list, detailed):
    align_output_list = align_calculate_RMSD(
        input_protein, chain_file_list, output_dir, chain_name_list)
    index = 0
    dict_unsorted = {}  # dictionary to store the chain_name:RMSD
    print(input_protein, " results of aligning")
    for align_output in align_output_list:
        print_align_output(align_output, chain_name_list[index], detailed)
        if (detailed):
            print("\n")
        build_dictionary(align_output, chain_name_list[index], dict_unsorted)
        index += 1
    print("\nMatching pdbs sorted by RMSD after aligning with ", input_protein)
    dict_sorted = get_RSMD_sorted(dict_unsorted)
    for k, v in dict_sorted.items():
        print(k, ":", v)
    return dict_unsorted, dict_sorted


detailed = False  # Change to True to see the detailed print output
cterminal_dict_unsorted, cterminal_dict_sorted = get_alignment_RMSDs(
    'c-terminal-rnf41', cTerminal_chain_file_list, cTerminal_output_dir, cTerminal_chain_name_list, detailed)
zf1_detailed = False  # Change to True to see the detailed print output
zf1_dict_unsorted, zf1_dict_sorted = get_alignment_RMSDs(
    'zf1-rnf41', zf1_chain_file_list, zf1_output_dir, zf1_chain_name_list, zf1_detailed)
zf2_detailed = False  # Change to True to see the detailed print output
zf2_dict_unsorted, zf2_dict_sorted = get_alignment_RMSDs(
    'zf2-rnf41', zf2_chain_file_list, zf2_output_dir, zf2_chain_name_list, zf2_detailed)


def pdb_CheckZincAtom_CheckResidues(base_structure_name, chain_name, chain_file):
    p = PDBParser()
    structure_name = base_structure_name + '-Aligned-' + chain_name
    structure = p.get_structure(structure_name, chain_file)
    has_zinc = False
    zinc_residues_list = []
    for residue in structure.get_residues():
        resname = residue.get_resname()
        if 'ZN' in resname:
            has_zinc = True
            res_id = resname + '-' + str(residue.get_full_id()[3][1])
            zinc_residues_list.append(res_id)
    if has_zinc:
        return has_zinc, zinc_residues_list
    cys_his_residues_list = []
    for residue in structure.get_residues():
        resname = residue.get_resname()
        if resname == 'CYS' or resname == 'HIS':
            res_id = resname + '-' + str(residue.get_full_id()[3][1])
            cys_his_residues_list.append(res_id)
    return has_zinc, cys_his_residues_list


def get_distances(chain_file, only_chain):
    parser = PDBParser()
    structure = parser.get_structure('TEST', chain_file)
    model = structure[0]
    chain = model[only_chain]
    distances_dict = {}
    match = False
    for residue1 in chain:
        for residue2 in chain:
            if residue1 != residue2:
                if 'HIS' in residue1.get_resname():
                    if 'HIS' in residue2.get_resname():
                        distance = residue1['CA'] - residue2['CA']
                        if distance <= 5:
                            match = True
                            res1 = residue1.get_resname() + '-' + \
                                str(residue1.get_full_id()[3][1])
                            res2 = residue2.get_resname() + '-' + \
                                str(residue2.get_full_id()[3][1])
                            distance_key = res1 + ' to ' + res2
                            distances_dict.update({distance_key: distance})
                    if 'CYS' in residue2.get_resname():
                        distance = residue1['CA'] - residue2['CA']
                        if distance <= 5:
                            match = True
                            res1 = residue1.get_resname() + '-' + \
                                str(residue1.get_full_id()[3][1])
                            res2 = residue2.get_resname() + '-' + \
                                str(residue2.get_full_id()[3][1])
                            distance_key = res1 + ' to ' + res2
                            distances_dict.update({distance_key: distance})
    return distances_dict, match


def check_zinc_cis_hys_distances(input_protein, chain_file_list, chain_name_list, dict_sorted):
    chain_name_index = 0
    zinc_dict = {}
    non_zinc_dict = {}
    distances_dict = {}
    for chain_file in chain_file_list:
        chain_name = chain_name_list[chain_name_index]
        has_zinc, residue_list = pdb_CheckZincAtom_CheckResidues(
            input_protein, chain_name, chain_file)
        if has_zinc:
            zinc_dict.update({chain_name: residue_list})
        else:
            only_chain = chain_name.split('-')
            distances_dict, match = get_distances(chain_file, only_chain[1])
            if match == True:
                if len(distances_dict) > 1:
                    dict_keys = list(distances_dict.keys())
                    for i in range(len(dict_keys)):
                        j = i + 1
                        if j < len(dict_keys):
                            key1_split = dict_keys[i].split(' to ')[0]
                            key2_split = dict_keys[j].split(' to ')[0]
                            if key1_split == key2_split:
                                non_zinc_dict.update(
                                    {chain_name: [residue_list, 'ZF_POSSIBLE', distances_dict]})
                                break
                else:
                    non_zinc_dict.update(
                        {chain_name: [residue_list, 'ZF_NOT_POSSIBLE', distances_dict]})

            else:
                non_zinc_dict.update(
                    {chain_name: [residue_list, 'ZF_NOT_POSSIBLE', distances_dict]})
        chain_name_index += 1
    dict_rsmd_residues = {}
    for key, value in dict_sorted.items():
        if key in zinc_dict:
            dict_rsmd_residues.update(
                {key: [value, 'ZINC_PRESENT', zinc_dict[key]]})
        elif key in non_zinc_dict:
            dict_rsmd_residues.update(
                {key: [dict_sorted[key], 'ZINC_ABSENT', non_zinc_dict[key]]})
        else:
            print("Error! Should not come here.")
    return dict_rsmd_residues

cterminal_dict_rmsd_residues = check_zinc_cis_hys_distances(
    'cterminal-rnf41', cTerminal_chain_file_list, cTerminal_chain_name_list, cterminal_dict_sorted)
print("CTERMINAL to matching proteins")
for k,v in cterminal_dict_rmsd_residues.items():
    print(k,":", v)

match = False
for k,v in cterminal_dict_rmsd_residues.items():
    if v[1] == 'ZINC_ABSENT' and v[2][1] == 'ZF_POSSIBLE':
        distances_dict = v[2][2] 
        if (len(distances_dict) > 1):
            print(k,":", v)
            match = True
    elif v[1] == 'ZINC_PRESENT':
        print(k,":ZINC ATOM PRESENT")
if match:
    print("Zinc Fingers possible. Check the above protein(s)")
else:
    print("ZINC FINGERS NOT A POSSIBILITY in REST OF PROTEINS")

zf1_dict_rmsd_residues_new = check_zinc_cis_hys_distances(
    'zf1-rnf41', zf1_chain_file_list, zf1_chain_name_list, zf1_dict_sorted)
print("ZF1 to matching proteins")
for k,v in zf1_dict_rmsd_residues_new.items():
    print(k,":", v)

match = False
for k,v in zf1_dict_rmsd_residues_new.items():
    if v[1] == 'ZINC_ABSENT' and v[2][1] == 'ZF_POSSIBLE':
        distances_dict = v[2][2] 
        if (len(distances_dict) > 1):
            print(k,":", v)
            match = True
    elif v[1] == 'ZINC_PRESENT':
        print(k,":ZINC ATOM PRESENT")
if match:
    print("Zinc Fingers possible. Check the above protein(s)")
else:
    print("ZINC FINGERS NOT A POSSIBILITY in REST OF PROTEINS")

zf2_dict_rmsd_residues = check_zinc_cis_hys_distances(
    'zf2-rnf41', zf2_chain_file_list, zf2_chain_name_list, zf2_dict_sorted)
%%capture cap --no-stderr
for k,v in zf2_dict_rmsd_residues.items():
    print(k,":", v)
with open('zf2-output.csv', 'w') as f:
    f.write(cap.stdout)
%%capture cap --no-stderr
for k,v in zf1_dict_rmsd_residues.items():
    print(k,":", v)
with open('zf1-output.csv', 'w') as f:
    f.write(cap.stdout)
%%capture cap --no-stderr
for k,v in cterminal_dict_rmsd_residues.items():
    print(k,":", v)
with open('cterminal-output.csv', 'w') as f:
    f.write(cap.stdout)
match = False
for k,v in zf2_dict_rmsd_residues.items():
    if v[1] == 'ZINC_ABSENT' and v[2][1] == 'ZF_POSSIBLE':
        distances_dict = v[2][2] 
        if (len(distances_dict) > 1): 
            print(k,":", v)
            match = True
    elif v[1] == 'ZINC_PRESENT':
        print(k,":matched")
if match:
    print("Zinc Fingers possible. Check the above protein(s)")
else:
    print("ZINC FINGERS NOT A POSSIBILITY in REST OF PROTEINS")

def common_data(nsp15File, rnf41File, list1, list2): 
    result = False
    for x in list1: 
        for y in list2: 
            p1 = x.split('-')[0]
            p2 = y.split('-')[0]
            if p1 == p2: 
                result = True
                print("NSP15:", nsp15File, " NSP15-matched-Protein-chain:", x, " NSP15-matched-Protein:", p1, "\n"
                      "RNF41-Protein:", rnf41File," RNF41-matched-Protein-chain:", y, " RNF41-matched-Protein:", p2)  
            if x == y: 
                result = True
                print("Chain match:", x, y)              
    if not result:
        print("\n\nNo matches for ", "NSP15:", nsp15File, "RNF41-Protein:", rnf41File, "\n")
    return result

def common_proteins(nsp15File): 
    cterminalFile = './CTerminalData.csv'
    zf1File = './zf1.csv'
    zf2File = './zf2.csv'
    nsp15DF = pd.read_csv(nsp15File)
    nsp15List = nsp15DF["Chains"].tolist()
    cTerminalDF = pd.read_csv(cterminalFile)
    cTerminalList = cTerminalDF["Chains"].tolist()
    common_data(nsp15File, cterminalFile, nsp15List, cTerminalList)
    zf1DF = pd.read_csv(zf1File)
    zf1List = zf1DF["Chains"].tolist()
    common_data(nsp15File, zf1File, nsp15List, zf1List)
    zf2DF = pd.read_csv('./zf2.csv')
    zf2List = zf2DF["Chains"].tolist()
    common_data(nsp15File, zf2File, nsp15List, zf2List)

result = common_proteins('./2gwf-usp8-chains.csv')
another_result = common_proteins('./6vww-nsp15-chains.csv')

nsp15DF = pd.read_csv('./6vww-nsp15-chains.csv') #nsp15File
nsp15List = nsp15DF["Chains"].tolist()

usp8DF = pd.read_csv('2gwf-usp8-chains.csv') #usp8 file
usp8List = usp8DF["Chains"].tolist()

third_result = common_data('./6vww-nsp15-chains.csv','./2gwf-usp8-chains.csv' , nsp15List, usp8List)
cmd.load('./6vww.pdb')
cmd.load('./zf1-RNF41-Dali-Output-active-sites-pdb/5vgc-A.pdb')
output = cmd.align('6vww', '5vgc-A')
print(output)
