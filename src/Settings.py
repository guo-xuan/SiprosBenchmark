'''
Created on Jun 27, 2017

@author: xgo
'''

import sys

label_train_str = 'Rev_' # 'Rev_'
label_test_str = 'TestRev_' # 'Shu_'
label_reserve_str = 'Rev_2_'

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

# remove protein not forward and testing
def clean_protein(protein_list):
    new_list = []
    for protein_str in protein_list:
        if not (protein_str.startswith(label_train_str)):
            if protein_str not in new_list:
                new_list.append(protein_str)
    return new_list


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