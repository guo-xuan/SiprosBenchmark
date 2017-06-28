'''
Created on Jun 28, 2017

@author: xgo
'''

import re, sys

import Settings

pattern = re.compile('[\W_]+')

def get_pep_from_pin(input_file_str, output_file_str):
    
    # input_file = '/media/xgo/Seagate/Proteomics/Experiments/BenchmarkRev/Comet/D10/pin/all_re_index.pin'
    # out_file = '/media/xgo/Seagate/Proteomics/Experiments/BenchmarkRev/Comet/D10/pin/pep_before_filtering.pin'
    fw = open(output_file_str, 'w')
    min_int = 10000
    max_int = 0
    pep_set = set()
    with open(input_file_str, 'r') as fr:
        header_str = fr.readline()
        header_list = header_str.strip().split('\t')
        pep_idx = header_list.index("peptide")
        for line_str in fr:
            split_list = line_str.strip().split('\t')
            pep_str = split_list[pep_idx]
            pep_str = pep_str[2:-2]
            pep_str = pattern.sub('', pep_str)
            if pep_str in pep_set:
                continue
            else:
                pep_set.add(pep_str)
            fw.write(pep_str)
            fw.write('\n')
            if len(pep_str) < min_int:
                min_int = len(pep_str)
            if len(pep_str) > max_int:
                max_int = len(pep_str)
            
    fw.close()
    print("max\t%d" % max_int)
    print("min\t%d" % min_int)
    print("Done.")
    
def complete_proteins_and_clean(input_file_str, pep_sip_file_str, output_file_str):
    
    fw = open(output_file_str, 'w')
    pep_pro_dict = {}
    with open(pep_sip_file_str, 'r') as fr:
        for line_str in fr:
            split_list = line_str.strip().split('\t')
            if len(split_list) < 2: # peptide too short or too long not in the sip file
                continue
            pro_list = split_list[1:]
            if split_list[0] in pep_pro_dict:
                l = pep_pro_dict[split_list[0]]
                for pro in pro_list:
                    if pro not in l:
                        l.append(pro)
            else:
                l = []
                pep_pro_dict[split_list[0]] = l
                for pro in pro_list:
                    if pro not in l:
                        l.append(pro)
    with open(input_file_str, 'r') as fr:
        header_str = fr.readline()
        fw.write(header_str)
        header_list = header_str.strip().split('\t')
        pep_idx = header_list.index("peptide")
        for line_str in fr:
            split_list = line_str.strip().split('\t')
            pep_str = split_list[pep_idx]
            pep_str = pep_str[2:-2]
            pep_str = '[' + pattern.sub('', pep_str) + ']'
            if len(pep_str) < 9 or len(pep_str) > 62:
                continue
            if pep_str in pep_pro_dict:
                pro_list = pep_pro_dict[pep_str]
            else:
                pro_list = split_list[(pep_idx+1):]
            if len(pro_list) == 0:
                print('error' + pep_str)
            protein_type_int = Settings.get_protein_type(pro_list)
            if protein_type_int == Settings.LabelRevTrain or protein_type_int == Settings.LabelRevReserve:
                split_list[1] = '-1'
            else:
                split_list[1] = '1'
            fw.write('\t'.join(split_list[0:(pep_idx+1)]))
            fw.write('\t')
            fw.write('\t'.join(pro_list))
            fw.write('\n')
    fw.close()
    print('complete_proteins_and_clean is done.')
    
def main(args=None):
    '''
    # extract peptides from pin files
    input_file_str = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/Ecoli/Results/Comet/MS2_marine/pin/all_new.pin"
    output_file_str = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/Ecoli/Results/Comet/MS2_marine/pin/all_new_pep.txt"
    get_pep_from_pin(input_file_str, output_file_str)
    '''
    
    # complete the proteins and remove too short or too long peptides, also correct lables
    input_file_str = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/Ecoli/Results/Comet/MS2_marine/pin/all_new.pin"
    pep_sip_file_str = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/Ecoli/Results/Comet/MS2_marine/pin/all_new_pep.sip"
    output_file_str = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/Ecoli/Results/Comet/MS2_marine/pin/all_corrected.pin"
    complete_proteins_and_clean(input_file_str, pep_sip_file_str, output_file_str)
    
    
if __name__ == '__main__':
    sys.exit(main())
    