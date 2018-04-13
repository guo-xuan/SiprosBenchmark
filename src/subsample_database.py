'''
Created on Jul 19, 2017

@author: xgo
'''

import random
import sys

def subsample(filename_str: str, output_str: str, subsample_rate: float):
    seq_str = ''
    id_str = ''
    with open(output_str, 'w') as fw:
        with open(filename_str, 'r') as fr:
            for line_str in fr:
                if line_str.startswith('>'):
                    if seq_str != '':
                        if random.random() < subsample_rate:
                            fw.write(id_str)
                            fw.write(seq_str)
                        id_str = ''
                        seq_str = ''
                    id_str = line_str
                else:
                    seq_str += line_str
    print('subsample is done.')
                    
def main():
    input_file_str = '/media/xgo/Seagate/Proteomics/Data/Angelo/Database/13_combo.fasta'
    output_file_str = '/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/Ecoli/Data/DB/13_combo_subsample_0.1.fasta'
    subsample_rate = 0.1
    subsample(input_file_str, output_file_str, subsample_rate)

if __name__ == '__main__':
    sys.exit(main())