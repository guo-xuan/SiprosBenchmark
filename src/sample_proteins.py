'''
Created on Jul 12, 2017

@author: xgo
'''

import sys
import random

def sample_protein_database(filename_str: str, sample_rate_float: float, output_file_str: str) -> None:
    keep_sequence_bool = False
    with open(output_file_str, 'w') as fw:
        with open(filename_str, 'r') as fr:
            for line_str in fr:
                if line_str.startswith('>'):
                    if random.random() < sample_rate_float:
                        keep_sequence_bool = True
                        fw.write(line_str)
                    else:
                        keep_sequence_bool = False
                else:
                    if keep_sequence_bool == True:
                        fw.write(line_str)
    print('sample_protein_database is done.')
    
def main():
    
    filename_str = '/media/xgo/Seagate/Proteomics/Data/Ecoli/Ecoli_K12_MG1655.fasta'
    sample_rate_float = 0.5
    output_file_str = '/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/Ecoli/Data/DB/Ecoli_K12_MG1655_sample_0.5.fasta'
    
    sample_protein_database(filename_str, sample_rate_float, output_file_str)

if __name__ == '__main__':
    sys.exit(main())