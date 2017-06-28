#!/usr/bin/python

## Import standard modules
from __future__ import division
import sys, getopt,  re
from collections import defaultdict
import csv 

## Import Sipros package modules
import parseconfig
import Settings

## increase CSV field size limit
csv.field_size_limit(1000000000)

## Version control
def get_version():
    return "4.0.1 (Alpha)"
"""
1. for Sipros4.0 (added two columns in sip file)
   new program should handle both types of sip files
"""



## Import classes

## Exit system with error message
die = Settings.die

## Check file exist
check_file_exist = Settings.check_file_exist

## Get file(s) list in working dir with specific file extension
get_file_list_with_ext = Settings.get_file_list_with_ext

## Get base_out filename
get_base_out = Settings.get_base_out

## Class for PepOutFields object
PepOutFields = Settings.PepOutFields

## Class for PsmOutFields object
PsmOutFields = Settings.PsmOutFields
Psm4OutFields = Settings.Psm4OutFields

## list_to_string
list_to_string = Settings.list_to_string

## Division error handling
divide = Settings.divide

## PrettyFloat (%0.5f)
PrettyFloat = Settings.PrettyFloat

## Get item list
check_sub_list = Settings.check_sub_list

## Get item list
get_item_list = Settings.get_item_list


## Help message
help_message = '''

Usage:
    python sipros_peptides_assembling.py [options]

Inputs:
    sipros_peptides_filtering output: [*.pep.txt] and [*.psm.txt] files
        (multiple peptide filtering results can be processed)
        (search automatically in current directory)
    sipros config file    (search automatically in current directory)

Options:
    -h/--help
    -v/--version
    -w/--working-dir ./path/    # Directory path containing SIPROS output files
                                # (default = current directory) 
    -c/--config-file SiprosConfig.cfg    # SIPROS config file 
                                         # (default = SiprosConfig.cfg) 

Outputs:
    BaseFilename.pro.txt
    BaseFilename.pro2pep.txt
    BaseFilenameepro2psm.txt
        - where BaseFilename = common prefix of inputs
        - if len(BaseFilename) < 5, then BaseFilename = Sipros_searches
'''

## Glboal variables
pep_file_ext = '.pep.txt'
psm_file_ext = '.psm.txt'

pep_iden_str = '[Peptide_Identification]'
fasta_database_str = 'FASTA_Database'
pro_iden_str = '[Protein_Identification]'
decoy_prefix_str = 'Decoy_Prefix'
min_peptide_per_protein_str = 'Min_Peptide_Per_Protein'
min_unique_peptide_per_protein_str = 'Min_Unique_Peptide_Per_Protein'
remove_decoy_identification_str = 'Remove_Decoy_Identification'

sipros4_psmout_column_length = 20
sipros4_input = None

# defaul value
decoy_prefix = 'Rev_'
min_peptide_per_protein = 2
min_unique_peptide_per_protein = 1
remove_decoy_identification = 'No'


## Parse options
def parse_options(argv):

    opts, _args = getopt.getopt(argv[1:], "hvVw:c:",
                                    ["help",
                                     "version",
                                     "working-dir",
                                     "config-file"])

    # Default working dir and config file
    working_dir = "./"
    config_file = "SiprosConfig.cfg"

    # Basic options
    for option, value in opts:
        if option in ("-v", "-V", "--version"):
            print("sipros_peptides_assembling.py V%s" % (get_version()))
            sys.exit(0)
        if option in ("-w", "--working-dir"):
            working_dir = value
            if working_dir[-1] != '/':
                working_dir = working_dir + '/'
        if option in ("-c", "--config-file"):
            config_file = value

    # only -w is provided
    if working_dir != "./" and config_file == "SiprosConfig.cfg":
        config_file = working_dir + config_file

    return (working_dir, config_file)


## Parse config file
def parse_config(config_filename):

    # Save config values to dictionary
    config_dict = {}    # initialize dictionay

    # Call Yinfeng's parseconfig.py module
    check_file_exist(config_filename)
    all_config_dict = parseconfig.parseConfigKeyValues(config_filename)

    # only save protein_identification config info to config_dict
    config_dict[decoy_prefix_str] = decoy_prefix
    config_dict[min_peptide_per_protein_str] = min_peptide_per_protein
    config_dict[min_unique_peptide_per_protein_str] = min_unique_peptide_per_protein
    config_dict[remove_decoy_identification_str] = remove_decoy_identification
    for key, value in all_config_dict.items():
        if key == (pep_iden_str + fasta_database_str):
            config_dict[fasta_database_str] = value
        elif key == (pro_iden_str + decoy_prefix_str):
            config_dict[decoy_prefix_str] = value
        elif key == (pro_iden_str + min_peptide_per_protein_str):
            config_dict[min_peptide_per_protein_str] = value
        elif key == (pro_iden_str + min_unique_peptide_per_protein_str):
            config_dict[min_unique_peptide_per_protein_str] = value
        elif key == (pro_iden_str + remove_decoy_identification_str):
            config_dict[remove_decoy_identification_str] = value
        else:
            continue

    # return config dictionary
    return config_dict


## Read fasta file and save the description
def read_fasta_file(working_dir, config_dict):

    # get fasta filename from config file
    fasta_filename = config_dict[fasta_database_str].strip()
    fasta_filename_only = fasta_filename.split("/")[-1]

    # get working_dir
    if working_dir[-1] != '/':
        working_dir = working_dir + '/'

    fasta_filename_dir = working_dir + fasta_filename_only

    # check file exist and open the file
    try:
        with open(fasta_filename) as _f: pass
        fasta_file = open(fasta_filename, 'r')
    except:
        try:
            with open(fasta_filename_only) as _f: pass
            fasta_file = open(fasta_filename_only, 'r')
        except:
            try:
                with open(fasta_filename_dir) as _f: pass
                fasta_file = open(fasta_filename_dir, 'r')
            except:
                print >> sys.stderr, '\nCannot open', fasta_filename
                print >> sys.stderr, 'Check your config file!'
                die("Program exit!")

    # save the fasta ID and description to the fasta_ID_dict
    fasta_ID_dict = {}    # initialize dictionary

    # FASTA ID is space delimited file
    fasta_ID_del = ' '

    # read lines
    for line in fasta_file.readlines():
        line = line.strip()
        # check line start with '>'
        if line.startswith('>'):
            # protein ID is the first word without '>'
            protein_ID = line.split(fasta_ID_del)[0][1:]
            # protein description is the whole line without '>'
            protein_desc = line[1:]
            # replace "#" to "$"
            protein_desc = re.sub('#', '$', protein_desc)
            # replace tab to " "
            protein_desc = re.sub('\t', ' ', protein_desc)

            # save the fasta ID and description to the fasta_ID_dict
            fasta_ID_dict[protein_ID] = protein_desc    # initialize dictionary

    return fasta_ID_dict


## check pep and psm files pair set, and save run#
def get_run_num(pep_file_list, psm_file_list):

    # dictionary of Run# for each pep_file
    run_num_dict = {}
    psm_run_num_dict = {}

    # read multiple pep files
    for pep_file_idx, pep_file in enumerate(pep_file_list):
    
        # check file exist and open
        check_file_exist(pep_file)

        # If base common prefix ends with '.pep.txt', then remove '.pep.txt'
        base_pep_file = pep_file.replace(pep_file_ext, "")

        # make a psm filename using pep filename
        psm_file = base_pep_file + psm_file_ext

        # check psm file exist
        check_file_exist(psm_file)

        # for Run#
        run_tag = 'Run' + str(pep_file_idx + 1)

        # run_num_dict
        run_num_dict[pep_file] = run_tag
        psm_run_num_dict[psm_file] = run_tag

    return (run_num_dict, psm_run_num_dict)

def get_pep_score_targetmatch(run_num_dict):
    pep_score_targetmatch_list = []
    # key = pep_file, val = run_num , sorted by Run# index
    for pep_file, _run_num in sorted(run_num_dict.items(), key=lambda x: x[1][-1]):
        # read line with csv
        fr = open(pep_file, 'r')
        
        # skip header
        fr.readline()

        # get data
        for pep_line in fr:
            
            if pep_line.startswith('#') or pep_line.startswith("IdentifiedPeptide"):
                continue
            
            # adapt class PepOutFields
            pep_list = pep_line.split('\t')
            pep_obj = PepOutFields._make(pep_list)
            l = [float(pep_obj.BestScore), pep_obj.TargetMatch]
            pep_score_targetmatch_list.append(l)
            
        fr.close()
    
    pep_score_targetmatch_sorted_list = sorted(pep_score_targetmatch_list, key = lambda x: x[0], reverse = True)
    
    return pep_score_targetmatch_sorted_list

# # FDR calculator
def FDR_calculator(FP, TP):
    FDR_numerator = float(FP)
    FDR_denominator = float(TP)
    FDR_accept = True

    if  FDR_denominator == 0:
        FDR_value = 1.0
        FDR_accept = False
    else:
        FDR_value = divide(FDR_numerator, FDR_denominator)
        FDR_accept = True

    return (FDR_accept, float(FDR_value))

def get_pep_score_cutoff(pep_score_targetmatch_sorted_list, fdr):
    fwr_int = 0
    rev_int = 0
    best_fwr = 0
    best_shu = 0
    score_cutoff = 0
    for one in pep_score_targetmatch_sorted_list:
        if one[1] == 'T':
            fwr_int += 1
        else:
            rev_int += 1
        
        (FDR_accept, FDR_value) = FDR_calculator(rev_int, fwr_int)
        if (FDR_accept is True) and (FDR_value <= fdr) and ((fwr_int + rev_int) > (best_fwr + best_shu)) :
            best_fwr = fwr_int
            best_shu = rev_int
            score_cutoff = one[0]
    
    return score_cutoff

# Read and load pep and psm files
def read_run_files(run_num_dict, pep_score_cutoff=None):

    # save the pep file data to the defaultdict
    pep_data_dict = defaultdict(list)
    # save the psm file data to the defaultdict
    psm_data_dict = defaultdict(list)
    # save peptides list (value) to unique protein (key)
    pro_pep_dict = defaultdict(list)
    # save protein list (value) to unique peptide (key)
    pep_pro_dict = defaultdict(list)

    # key = pep_file, val = run_num , sorted by Run# index
    for pep_file, run_num in sorted(run_num_dict.items(), key=lambda x: x[1][-1]):

        # read line with csv
        pep_reader = open(pep_file, 'r')
        # skip header
        pep_reader.readline()

        # get data
        for pep_line in pep_reader:
            
            if pep_line.startswith('#') or pep_line.startswith("IdentifiedPeptide"):
                continue
            
            pep_line = pep_line.split('\t')
            
            # adapt class PepOutFields
            pep_obj = PepOutFields._make(pep_line)
            if pep_score_cutoff is not None and float(pep_obj.BestScore) < pep_score_cutoff:
                continue
            identified_peptide = pep_obj.IdentifiedPeptide
            identified_peptide = identified_peptide.strip()
         
            # new unique ID for pep_data_dict is identified_peptide + run_num
            pep_run_id = identified_peptide + "_" + run_num

            pep_data_dict[pep_run_id].append(pep_line)

            # get protein item list
            protein_names = pep_obj.ProteinNames
            pro_item_list = get_item_list(protein_names.strip())

            # for loop of pro_item_list
            for pro_item in pro_item_list:

                pro_item = pro_item.strip()
                
                # save IdentifiedPeptide list to unique protein dict
                if identified_peptide not in pro_pep_dict[pro_item]:
                    pro_pep_dict[pro_item].append(identified_peptide)

                # save protein list to unique peptide dict
                if pro_item not in pep_pro_dict[identified_peptide]:
                    pep_pro_dict[identified_peptide].append(pro_item)
                    
        pep_reader.close()

    # try to find indistinguishable set
    indistin_pro_dict = defaultdict(list)
    for pro_key, pep_list in pro_pep_dict.items():
        sorted_pep_list = sorted(set(pep_list))
        sorted_pep_list_join = '_'.join(sorted_pep_list)
        indistin_pro_dict[sorted_pep_list_join].append(pro_key)

    # indistin_key = str(sorted(pep_list)), indistin_value=pro_key list
    for _indistin_key, indistin_value in indistin_pro_dict.items():
        # if proteins have a same set of peptides

        if len(indistin_value) > 1:
            # get new protein name
            new_pro_key = list_to_string(sorted(indistin_value))
            new_pro_val = pro_pep_dict[indistin_value[0]]
            # append with new key and value to pro_pep_dict
            # delete provious pep_pro_dict key,value, and append new one
            for new_pro_val_one in new_pro_val:
                pro_pep_dict[new_pro_key].append(new_pro_val_one)
                pep_pro_dict[new_pro_val_one].append(new_pro_key)

                # delete pep_pro_dict indistin protein list of for the peptide
                org_pep_pro_dict_list = pep_pro_dict[new_pro_val_one]
                for indistin_value_one in indistin_value:
                    org_pep_pro_dict_list.remove(indistin_value_one)

                if len(org_pep_pro_dict_list) == 0:
                    del pep_pro_dict[new_pro_val_one]
                else:
                    pep_pro_dict[new_pro_val_one] = org_pep_pro_dict_list

            # delete previous indistinguishable proteins in pro_pep_dict
            for indistin_pro in indistin_value:
                del pro_pep_dict[indistin_pro]

    # key = pep_file, val = run_num , sorted by Run# index
    for pep_file, run_num in sorted(run_num_dict.items(), key=lambda x: x[1][-1]):

        # If base common prefix ends with '.pep.txt', then remove '.pep.txt'
        base_pep_file = pep_file.replace(pep_file_ext, "")

        # make a psm filename using pep filename
        psm_file = base_pep_file + psm_file_ext

        # check psm file exist
        check_file_exist(psm_file)

        # read line with csv
        psm_reader = open(psm_file, 'r')
        # skip header
        psm_reader.readline()

        # get data
        for psm_line in psm_reader:
            
            if psm_line.startswith('#') or psm_line.startswith("Filename"):
                continue
            
            psm_line = psm_line.split('\t')
            
            # check sipros3 or sipros4 input
            psm_obj = None 
            global sipros4_input
            if len(psm_line) == sipros4_psmout_column_length:
                psm_obj = Psm4OutFields._make(psm_line)
                sipros4_input = True
            else:
                psm_obj = PsmOutFields._make(psm_line)
                sipros4_input = False
            
            if pep_score_cutoff is not None and float(psm_obj.Score) < pep_score_cutoff:
                continue
            
            # get protein item list
            protein_names = psm_obj.ProteinNames
            pro_item_list = get_item_list(protein_names.strip())

            # save key=protein_id, val=line
            for pro_item_one in pro_item_list:

                # new unique ID for psm_data_dict is protein_name + run_num
                pro_item_id = pro_item_one + "_" + run_num

                # save to dictionary
                psm_data_dict[pro_item_id].append(psm_line)
        
        psm_reader.close()
        
    return (pep_data_dict, psm_data_dict, pro_pep_dict, pep_pro_dict)


## Greedy algorithm to extract proteins
def greedy_alg(config_dict, pro_pep_dict, pep_pro_dict):

    # get config value
    min_peptide_per_protein = int(config_dict[min_peptide_per_protein_str])
    min_unique_peptide_per_protein = int(config_dict[min_unique_peptide_per_protein_str])

    # return config dictionary

    # First, extract proteins that have >= min_peptide_per_proteins 
    # and at least min_unique_peptide_per_protein of those is unique

    # save the extracted proteins to the list
    pro_greedy_list = []

    # copy pro_pep_dict to pro_pep_dict_red for reduction in greedy steps
    pro_pep_dict_red = defaultdict(list)
    for pro_key, pep_list in pro_pep_dict.items():
        for pep_list_one in pep_list:
            pro_pep_dict_red[pro_key].append(pep_list_one)

    # copy pep_pro_dict to pep_pro_dict_red for reduction in greedy steps
    pep_pro_dict_red = defaultdict(list)
    for pep_key, pro_list in pep_pro_dict.items():
        for pro_list_one in pro_list:
            pep_pro_dict_red[pep_key].append(pro_list_one)

    # for loop of pro_pep_dict
    for pro_key, pep_list in pro_pep_dict.items():

        # if proteins that have >= min_peptide_per_protein
        if len(pep_list) >= min_peptide_per_protein:

            # if at least min_unique_peptide_per_protein is unique
            unique_pep_pro_num = 0
            for pep_list_one in pep_list:
                pep_pro_num = len(pep_pro_dict[pep_list_one])
                if pep_pro_num == 1:
                    unique_pep_pro_num += 1
            # unique peptides num should be >= min_unique_peptide_per_protein
            # if min_unique_peptide_per_protein = 0, then set 1
            if unique_pep_pro_num >= max(min_unique_peptide_per_protein,1):

                # append to the pro_greedy_list
                if pro_key not in pro_greedy_list:
                    pro_greedy_list.append(pro_key)

                # remove the protein
                try:
                    del pro_pep_dict_red[pro_key]
                except:
                    pass
            
                # remove all peptides that are covered by the protein
                for pep_list_item in pep_list:
                    # delete by key
                    try:
                        del pep_pro_dict_red[pep_list_item]
                    except:
                        pass

                    # get new peptide list for the proteins from pep_pro_dict
                    for pro_list_item in pep_pro_dict[pep_list_item]:
                        new_pep_list = pro_pep_dict[pro_list_item]

                        # if new_pep_list is sub list of pep_list, then remove
                        if check_sub_list(new_pep_list, pep_list):
                            try:
                                del pro_pep_dict_red[pro_list_item]
                            except:
                                pass

    # Second, iteratively extract a protein at a time that covers the most peptides
    if len(pro_pep_dict_red.keys()) > 0:
        # Run greedy iterations until it converges
        converge = False

        # Iterate greedy algorithm until it converges
        greedy_step = 0
        while (converge == False):

            greedy_step += 1

            # find a protein that covers the most peptides
            ppdr_idx = 0
            for key_ppdr, val_ppdr in pro_pep_dict_red.items():
                if ppdr_idx == 0:
                    max_key_ppdr = key_ppdr
                    max_len_val_ppdr = len(val_ppdr)
                else:
                    # get current one
                    cur_len_val_ppdr = len(val_ppdr)
                    cur_key_ppdr = key_ppdr
                    # get max one
                    if cur_len_val_ppdr > max_len_val_ppdr:
                        max_len_val_ppdr = cur_len_val_ppdr 
                        max_key_ppdr = cur_key_ppdr
                ppdr_idx += 1

            max_pro_one = max_key_ppdr

            # if proteins that have >= min_peptide_per_protein
            if len(pro_pep_dict_red[max_pro_one]) >= min_peptide_per_protein:

                # if at least min_unique_peptide_per_protein is unique
                unique_pep_pro_num = 0
                for pep_list_one in pro_pep_dict_red[max_pro_one]:
                    pep_pro_num = len(pep_pro_dict[pep_list_one])
                    if pep_pro_num == 1:
                        unique_pep_pro_num += 1
                # if at least min_unique_peptide_per_protein is unique
                if unique_pep_pro_num >= min_unique_peptide_per_protein:

                    # append the protein to the pro_greedy_list
                    pro_greedy_list.append(max_pro_one)

                    # loop for pep set
                    for pep_list_sub_one in pro_pep_dict[max_pro_one]:
                        # delete by key
                        try:
                            del pep_pro_dict_red[pep_list_sub_one]
                        except:
                            pass

                        # loop for pro set 
                        for pro_list_sub_one in pep_pro_dict[pep_list_sub_one]:
                            try:
                                del pro_pep_dict_red[pro_list_sub_one]
                            except:
                                pass

                    # remove the protein
                    try:
                        del pro_pep_dict_red[max_pro_one]
                    except:
                        pass
            
                    # remove all peptides that are covered by the protein
                    for pep_list_item in pro_pep_dict[max_pro_one]:
                        # delete by key
                        try:
                            del pep_pro_dict_red[pep_list_item]
                        except:
                            pass

                        # get new peptide list for the proteins from pep_pro_dict
                        for pro_list_item in pep_pro_dict[pep_list_item]:
                            new_pep_list = pro_pep_dict[pro_list_item]
        
                            # if new_pep_list is sub list of pep_list, then remove
                            if check_sub_list(new_pep_list, pep_list):
                                try:
                                    del pro_pep_dict_red[pro_list_item]
                                except:
                                    pass
                else: 
                    del pro_pep_dict_red[max_pro_one]
            # the max peptides number for protein < min_peptide_per_protein
            else:
                converge = True

            # if there is no protein, then converge
            if len(pro_pep_dict_red.keys()) == 0:
                converge = True

        # greedy algorithm done

    pro_greedy_list = sorted(pro_greedy_list)

    return pro_greedy_list


## Get protein description to handle multiple protein IDs
def get_protein_description(protein_ID, fasta_ID_dict):

    # initialize
    protein_description_list = []

    # if multiple IDs
    if (protein_ID.startswith('{')) and (protein_ID.endswith('}')):

        # multiple proteins exist
        protein_ID_list = get_item_list(protein_ID.strip())
        for protein_ID_one in protein_ID_list:
            
            # check protein ID exist
            if protein_ID_one in fasta_ID_dict:
                protein_description_one = fasta_ID_dict[protein_ID_one]
                protein_description_list.append(protein_description_one)
            else: 
                protein_description_one = "N/A"
                protein_description_list.append(protein_description_one)

    # single ProteinID
    else:
        # check protein ID exist
        if protein_ID in fasta_ID_dict:
            protein_description = fasta_ID_dict[protein_ID]
            protein_description_list.append(protein_description)
        else: 
            protein_description = "N/A"
            protein_description_list.append(protein_description)

    # convert list to string
    protein_description = list_to_string(protein_description_list)

    return protein_description


## check decoy match
def check_decoy_match(ProteinNames, decoy_prefix):

    # match type (correct or decoy) and strip
    match_type = ProteinNames
    match_type_list = []

    if match_type.startswith("{"):
        # if starts with {}, then delete parenthesis { }
        match_type = match_type[1:-1]
        # sometimes multiple proteins
        match_type_list = re.split(r"\s*[,]\s*", match_type.strip())
    else:
        match_type_list.append(match_type)
    # TF -> True or False(decoy match)
    TF = False

    # for loop of proteins
    for match_item in match_type_list:
        # if at least one of matches is True, then match is True
        if not match_item.startswith(decoy_prefix):
            TF = True 
            break

    return TF

## check decoy match
def check_target_match(ProteinNames, target_prefix):

    # match type (correct or decoy) and strip
    match_type = ProteinNames
    match_type_list = []

    if match_type.startswith("{"):
        # if starts with {}, then delete parenthesis { }
        match_type = match_type[1:-1]
        # sometimes multiple proteins
        match_type_list = re.split(r"\s*[,]\s*", match_type.strip())
    else:
        match_type_list.append(match_type)
    # TF -> True or False(decoy match)
    TP = False

    # for loop of proteins
    for match_item in match_type_list:
        # if at least one of matches is True, then match is True
        if match_item.startswith(target_prefix):
            TP = True 
            break

    return TP


## Report output files
def report_output(config_dict,
                  run_num_dict,
                  psm_run_num_dict,
                  pep_data_dict,
                  psm_data_dict,
                  pro_pep_dict,
                  pep_pro_dict,
                  pro_greedy_list,
                  fasta_ID_dict=None,
                  target_prefix = 'lcl'):

    # get config value
    min_peptide_per_protein = int(config_dict[min_peptide_per_protein_str])
    min_unique_peptide_per_protein = int(config_dict[min_unique_peptide_per_protein_str])

    # to get decoy_prefix
    decoy_prefix = Settings.label_test_str

    # total number of proteins
    total_proteins_before_filtering = len(pro_pep_dict.keys())
    decoy_proteins_before_filtering = 0
    for key, _val in pro_pep_dict.items():
        check_decoy_match_val = check_decoy_match(key, decoy_prefix)
        if check_decoy_match_val is False:
            decoy_proteins_before_filtering += 1
    _target_proteins_before_filtering = int(total_proteins_before_filtering) - int(decoy_proteins_before_filtering)

    # total number of identified proteins
    total_proteins_after_filtering = len(pro_greedy_list)
    decoy_proteins_after_filtering = 0
    true_target_proteins_after_filtering = 0
    for pro_one in pro_greedy_list:
        check_decoy_match_val = check_decoy_match(pro_one, decoy_prefix)
        if check_decoy_match_val is False:
            decoy_proteins_after_filtering += 1
        if target_prefix is not None:
            check_target_match_val = check_target_match(pro_one, target_prefix)
            if check_target_match_val is True:
                true_target_proteins_after_filtering += 1
    target_proteins_after_filtering = int(total_proteins_after_filtering) - int(decoy_proteins_after_filtering)

    # protein FDR
    protein_fdr = 0.0
    if float(total_proteins_after_filtering) != 0:
        protein_fdr = divide(float(decoy_proteins_after_filtering), float(total_proteins_after_filtering))
    protein_fdr = PrettyFloat(protein_fdr)
    
    print("{:,d} ({:.3f}%) ({:,d})".format(target_proteins_after_filtering, 3.0*100.0*float(decoy_proteins_after_filtering)/float(target_proteins_after_filtering), true_target_proteins_after_filtering))

    return ((3*100.0*float(decoy_proteins_after_filtering)/float(target_proteins_after_filtering)) , target_proteins_after_filtering, true_target_proteins_after_filtering)


## +------+
## | Main |
## +------+
def main(argv=None):
    # parse options
    # (working_dir, config_filename) = parse_options(argv)

    working_dir = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/Ecoli/Results/SiprosEnsemble/soil/SiprosEnsemblePercolator/fdr_0.02/"
    config_filename = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/Ecoli/Data/Configs/SiprosEnsemble/soil.cfg"

    config_dict = parse_config(config_filename)

    # Get .pep.txt output file(s) in working directory
    pep_file_list = get_file_list_with_ext(working_dir, pep_file_ext)
    # Get .psm.txt output file(s) in working directory
    psm_file_list = get_file_list_with_ext(working_dir, psm_file_ext)
    # check pep and psm files pair set, and save run#
    (run_num_dict, psm_run_num_dict) = get_run_num(pep_file_list, psm_file_list)
    # Get base_out for output
    base_out_default = 'Sipros_searches'
    _base_out = get_base_out(pep_file_list, base_out_default, working_dir)
        
    protein_num_list = []
    true_protein_num_list = []
            
    fdr_float = 0.0000
    pep_score_targetmatch_sorted_list = get_pep_score_targetmatch(run_num_dict)
    try_max = 10
    target_protein_fdr_list = [0.5, 1.0, 2.0]
    
    for target_protein_fdr in target_protein_fdr_list:
        increasing_rate = 0.003
        num_try = 0
        while True:
            fdr = fdr_float + increasing_rate
            score_cutoff = get_pep_score_cutoff(pep_score_targetmatch_sorted_list, fdr)
            (pep_data_dict, psm_data_dict, pro_pep_dict, pep_pro_dict) = read_run_files(run_num_dict, score_cutoff)
            (pro_greedy_list) = greedy_alg(config_dict, pro_pep_dict, pep_pro_dict)
            print(str(fdr))
            (protein_fdr, protein_num, true_protein_num) = report_output(config_dict,
                                            run_num_dict,
                                            psm_run_num_dict,
                                            pep_data_dict,
                                            psm_data_dict,
                                            pro_pep_dict,
                                            pep_pro_dict,
                                            pro_greedy_list)
            num_try += 1
            if abs(protein_fdr - target_protein_fdr) < 0.01 or num_try > try_max:
                print('Done')
                protein_num_list.append(protein_num)
                true_protein_num_list.append(true_protein_num)
                break
            elif protein_fdr > target_protein_fdr:
                increasing_rate /= 2.0
            else:
                fdr_float += increasing_rate
                
    ecoli_mode = True
    
    if ecoli_mode:
        for idx, protein_num in enumerate(protein_num_list):
            sys.stdout.write('{:,d} '.format(protein_num))
            sys.stdout.write('({:,d})\t'.format(true_protein_num_list[idx]))        
    else:
        for protein_num in protein_num_list:
            sys.stdout.write('{:,d}\t'.format(protein_num))
    
## If this program runs as standalone, then exit.
if __name__ == "__main__":
    sys.exit(main())

