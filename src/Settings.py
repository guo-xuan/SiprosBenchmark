'''
Created on Jun 27, 2017

@author: xgo
'''

import sys, os, re, math
from collections import namedtuple

label_train_str = 'Rev_'
label_test_str = 'TestRev_'
label_reserve_str = 'Rev_2_'
label_ecoli_str = 'lcl'

LabelEcoli = 0
LabelFwd = 1
LabelRevTrain = 2
LabelRevReserve = 3
LabelTest = 4

# # Decoy Reverse Forward protein
def get_protein_type(protein_split_list):
    for one_protein in protein_split_list:
        if not (one_protein.startswith(label_train_str) or one_protein.startswith(label_test_str)):
            return LabelFwd
    for one_protein in protein_split_list:
        if one_protein.startswith(label_test_str):
            return LabelTest
    if label_reserve_str != '':
        for one_protein in protein_split_list:
            if one_protein.startswith(label_reserve_str):
                return LabelRevReserve
    return LabelRevTrain

def get_protein_type_ecoli(protein_split_list):
    for one_protein in protein_split_list:
        if one_protein.startswith(label_ecoli_str):
            return LabelEcoli

    for one_protein in protein_split_list:
        if not (one_protein.startswith(label_train_str) or one_protein.startswith(label_test_str)):
            return LabelFwd
    
    for one_protein in protein_split_list:
        if one_protein.startswith(label_test_str):
            return LabelTest
        
    if label_reserve_str != '':
        for one_protein in protein_split_list:
            if one_protein.startswith(label_reserve_str):
                return LabelRevReserve
    return LabelRevTrain

# remove protein not forward and testing
def clean_protein(protein_list):
    new_list = []
    for protein_str in protein_list:
        if not (protein_str.startswith(label_train_str) or protein_str.startswith(label_reserve_str)):
            if protein_str not in new_list:
                new_list.append(protein_str)
    return new_list

# remove protein not forward and testing
def remove_duplicate_protein(protein_list):
    new_list = []
    for protein_str in protein_list:
        if protein_str not in new_list:
            new_list.append(protein_str)
    return new_list

def findOccurences(s, ch):
    return len([i for i, letter in enumerate(s) if letter == ch])

FDR_parameter = 1.0
# # FDR calculator
def FDR_calculator(FP, TP):
    FDR_numerator = float(FP) * float(FDR_parameter)
    FDR_denominator = float(TP)
    FDR_accept = True

    if  FDR_denominator == 0:
        FDR_value = 1.0
        FDR_accept = False
    else:
        FDR_value = divide(FDR_numerator, FDR_denominator)
        FDR_accept = True

    return (FDR_accept, float(FDR_value))

# # Exit system with error message
def die(msg=None):
    if msg is not None:
        print >> sys.stderr, msg
        sys.exit(1)

# # Division error handling
def divide(x, y):
    try:
        result = x / y
    except ZeroDivisionError as detail:
        sys.stderr.write('Handling run-time error:', detail)
        die('Program exit!')
    else:
        return result
    
# # Get file(s) list in working dir with specific file extension
def get_file_list_with_ext(working_dir, file_ext):

    # define sipros file extension
    file_list = []

    # working directory
    if os.path.exists(working_dir):
        for file_name in os.listdir(working_dir):

            # check the file extension
            if file_name.endswith(file_ext):
                file_path_name = working_dir + file_name
                file_list.append(file_path_name)

        if len(file_list) == 0:
            sys.stderr.write("\nCannot open %s file(s)." % (file_ext))
            die("Program exit!")
        file_list = sorted(file_list)

    else:
        sys.stderr.write("\nCannot open working directory", working_dir)
        die("Program exit!")

    return file_list

# # Get base output filename with input file list and base_out_default
def get_base_out(file_list, base_out_default, working_dir):

    # Get base output with common prefix
    base_out = os.path.commonprefix(file_list)
    base_out_filename = base_out.split('/')[-1]

    # If base common prefix ends with '.pep.txt', then remove '.pep.txt'
    base_out = base_out.replace(".pep.txt", "_")

    # If base_out file name is less than 5, then use default baseout
    if len(base_out_filename) < 5:
        base_out = working_dir + base_out_default

    # If base common prefix ends with '_' or '.', then remove
    base_out = base_out[:-1] if (base_out[-1] in ('_', '.')) else base_out

    return base_out

# # get item list from parenthesis string as {AA,BB}
def get_item_list(input_string):

    input_string = input_string[1:-1]
    item_list = re.split(r"\s*[,]\s*", input_string.strip())
    item_list_new = []
    for item_one in item_list:
        if item_one not in item_list_new:
            item_list_new.append(item_one)
    return item_list_new

# # check sub list
def check_sub_list(list_A, list_B):

    check_status = True

    for list_A_item in list_A:
        if list_A_item not in list_B:
            check_status = False
        else:
            continue

    return check_status

# # Class for pretty float
class PrettyFloat(float):
    def __repr__(self):
        return "%0.5f" % self
    
# # list_to_string
# # if single element, then just convert to string
# # if multiple elements, then bracket {A,B}
def list_to_string(input_list):

    if len(input_list) > 1:
        converted_str = '{' + ','.join(input_list) + '}'
    else:
        converted_str = ''.join(input_list)

    return converted_str

# # Class for PepOutFields object
class PepOutFields(namedtuple('PepOutFields',
        ['IdentifiedPeptide',  # 0
         'ParentCharge',  # 1
         'OriginalPeptide',  # 2
         'ProteinNames',  # 3
         'ProteinCount',  # 4
         'TargetMatch',  # 5
         'SpectralCount',  # 6 number of PSMs matched to this peptide
         'BestScore',  # 7 the highest score of those PSMs
         'PSMs',  # 8 a list of PSMs matched to this peptide. Use {Filename[ScanNumber],Filename[ScanNumber]} format
         'ScanType',  # 9 ScanType
         'SearchName'])):  # 10 SearchName

    def __init__(self):
        self.data = self
        
# # Class for PsmOutFields object
class PsmOutFields(namedtuple('PsmOutFields',
        ['Filename',
         'ScanNumber',
         'ParentCharge',
         'MeasuredParentMass',
         'CalculatedParentMass',
         'MassErrorDa',  # CalculatedParentMass - MeasuredParentMass
         'MassErrorPPM',  # MassErrorDa / CalculatedParentMass
         'ScanType',
         'SearchName',
         'ScoringFunction',
         'Score',
         'DeltaZ',  # The difference score between the rank 1 and 2
         'DeltaP',  # The difference score between isoform
         'IdentifiedPeptide',
         'OriginalPeptide',
         'ProteinNames',
         'ProteinCount',
         'TargetMatch'])):
    def __init__(self):
        self.data = self
        
# # Class for PsmOutFields object (for sipro4)
class Psm4OutFields(namedtuple('PsmOutFields',
        ['Filename',
         'ScanNumber',
         'ParentCharge',
         'MeasuredParentMass',
         'CalculatedParentMass',
         'MassErrorDa',  # CalculatedParentMass - MeasuredParentMass
         'MassErrorPPM',  # MassErrorDa / CalculatedParentMass
         'ScanType',
         'SearchName',
         'ScoringFunction',
         'Score',
         'DeltaZ',  # The difference score between the rank 1 and 2
         'DeltaP',  # The difference score between isoform
         'IdentifiedPeptide',
         'OriginalPeptide',
         'ProteinNames',
         'ProteinCount',
         'TargetMatch',
         'AveAtom',
         'StdAtom'])):
    def __init__(self):
        self.data = self
        
        
# # Check file exist
def check_file_exist(filename):

    try:
        with open(filename) as _f: pass
    except IOError as _e:
        print >> sys.stderr, '\nCannot open', filename
        die("Program exit!")

fNeutronMass = 1.00867108694132 # it is Neutron mass
def get_mass_diff(measured_mass, calculated_mass):
    fDiff = calculated_mass - measured_mass
    fTemp = fDiff
    fCeil = 0
    fFloor = 0
    if fDiff >= 0:
        fDiff = fTemp
        fCeil = math.ceil(fTemp)*fNeutronMass
        fFloor = math.floor(fTemp)*fNeutronMass
        if fFloor > fTemp:
            fFloor -= fNeutronMass
        if fCeil - fNeutronMass > fTemp:
            fCeil -= fNeutronMass
        if fTemp > fCeil - fTemp:
            fTemp = fCeil - fTemp
        if fDiff > fDiff - fFloor:
            fDiff = abs(fDiff - fFloor)
        if abs(fTemp) < abs(fDiff):
            fDiff = fTemp
            dM = -fTemp
        else:
            dM = fDiff
    else:
        fCeil = math.ceil(fDiff)*fNeutronMass
        if fCeil < fDiff:
            fCeil += fNeutronMass
        fFloor = math.floor(fDiff)*fNeutronMass
        if fFloor + fNeutronMass < fDiff:
            fFloor += fNeutronMass
        fDiff = fTemp
        if abs(fTemp) > fCeil - fTemp:
            fTemp = fCeil - fTemp
        if abs(fDiff) > fDiff - fFloor:
            fDiff = fDiff - fFloor
        fTemp = abs(fTemp)
        fDiff = abs(fDiff)
        if fTemp < fDiff:
            fDiff = fTemp
            dM = -fTemp
        else:
            dM = fDiff
            
    fMassDiff = fDiff
    return dM, fMassDiff


def get_num_missed_cleavage_sites(sIdentifiedSeq, sResiduesBeforeCleavage = "KR", sResiduesAfterCleavage = "ACDEFGHIJKLMNPQRSTVWY"):
    count_int = 0
    for i in range(len(sIdentifiedSeq) - 1):
        if sIdentifiedSeq[i] in sResiduesBeforeCleavage and sIdentifiedSeq[i + 1] in sResiduesAfterCleavage:
            count_int += 1
    return count_int