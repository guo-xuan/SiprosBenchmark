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

# TestRev should be postive label
# internal missed cleavage site should be less than 4
def correct_pin(input_file_str, output_file_str):
    with open(output_file_str, 'w') as fw:
        with open(input_file_str, 'r') as fr:
            fw.write(fr.readline()) # header 1
            fw.write(fr.readline()) # header 2
            for line_str in fr:
                split_list = line_str.split('\t')
                if int(split_list[27]) > 3: # missed cleavage sites
                    continue
                protein_list = split_list[29:]
                protein_type = get_protein_type(protein_list)
                if protein_type == LabelRevTrain or protein_type == LabelRevReserve:
                    split_list[1] = "-1"
                else:
                    split_list[1] = "1"
                fw.write('\t'.join(split_list))
    
    print("Correct Pin Done.")

def main():
    
    # correct the protein labels
    input_file_str = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/Ecoli/Results/MsGF/ecoli_new/all.pin"
    output_file_str = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/Ecoli/Results/MsGF/ecoli_new/all_corrected.pin"
    correct_pin(input_file_str, output_file_str)


if __name__ == '__main__':
    sys.exit(main())