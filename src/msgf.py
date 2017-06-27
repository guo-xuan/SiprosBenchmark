'''
Created on Jun 27, 2017

@author: xgo
'''

import sys

import Settings

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
                protein_type = Settings.get_protein_type(protein_list)
                if protein_type == Settings.LabelRevTrain or protein_type == Settings.LabelRevReserve:
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