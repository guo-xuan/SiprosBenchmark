'''
Created on Jul 6, 2017

@author: xgo
'''
import sys
from os import listdir
from os.path import isfile, join

## global variables
## training prefix
train_str = 'Rev_1_' 
## testing prefix
test_str = 'TestRev_' 
## reserved prefix
reserve_str = 'Rev_2_'
## ratio of testing decoy vs forward
Test_Fwd_Ratio = 1
## maximum precurer mass windows
mass_window_max_int = 0
## feature list
feature_selection_list = []
## fwd psm value
LabelFwd = 1
## training psm value
LabelTrain = 2
## testing psm value
LabelTest = 3
## reserved psm value
LabelReserve = 4
## sip mode training psm value
LabelSipTrainFwd = 1
## forward psms
num_forward_psms_before_filtering = 0
## protein database size
num_proteins = 0

def get_protein_type(protein_sequence, lProtein=None):
    """
    get the protein type
    if all reserved type, return LabelReserve
    if all testing type, return LabelTest
    if all training type, return LabelTrain
    otherwise, it is forward protein, return LabelFwd
    """
    sProteins = protein_sequence.replace('{', '')
    sProteins = sProteins.replace('}', '')
    asProteins = sProteins.split(',')
    if lProtein != None:
        del lProtein[:]
        for sProtein in asProteins:
            sProtein = sProtein.strip()
            if sProtein not in lProtein:
                lProtein.append(sProtein)
                
    protein_list_tmp_1 = []
    protein_list_tmp_2 = []
    
    reserve_type = True
    if reserve_str != '':
        for sProtein in asProteins:
            if not sProtein.startswith(reserve_str):
                protein_list_tmp_1.append(sProtein)
                reserve_type = False
        if reserve_type:
            return LabelReserve
    
    training_type = True
    if train_str != '':
        for sProtein in protein_list_tmp_1:
            if not sProtein.startswith(train_str):
                protein_list_tmp_2.append(sProtein)
                training_type = False
        if training_type:
            return LabelTrain
    
    testing_type = True
    if test_str != '':
        for sProtein in protein_list_tmp_2:
            if not sProtein.startswith(test_str):
                testing_type = False
        if testing_type:
            return LabelTest
    
    return LabelFwd
                        
def main():
    input_file_str = '/media/naux/naux_data/proteomics/experiments/sipros_ensemble/ecoli_samples/sipros_ensemble/soil/EColi_Try_HCD_DE10ppm_CS_1000_NCE30_180_Run.tab'
    num_forward = 0.0
    num_total = 0.0
    with open(input_file_str, 'r') as fr:
        line_str = fr.readline()
        split_list = line_str.strip().split('\t')
        protein_idx = split_list.index('ProteinNames')
        for line_str in fr:
            split_list = line_str.strip().split('\t')
            protein_names = split_list[protein_idx]
            if get_protein_type(protein_names) == LabelFwd:
                num_forward += 1
            num_total += 1
    
    print(str(num_forward))
    print(str(num_total))
    print(str(num_forward/num_total))
    
if __name__ == '__main__':
    sys.exit(main())