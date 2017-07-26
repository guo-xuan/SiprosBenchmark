'''
Created on Jun 28, 2017

@author: xgo
'''

import re, sys

import Settings
from sipros import protein_type
from os import listdir
from os.path import isfile, join

pattern = re.compile('[\W_]+')

def get_pep_from_pep_txt(input_file_str, output_file_str):
    
    fw = open(output_file_str, 'w')
    min_int = 10000
    max_int = 0
    pep_set = set()
    pep_idx = 2
    with open(input_file_str, 'r') as fr:
        _header_str = fr.readline()
        for line_str in fr:
            split_list = line_str.strip().split('\t')
            pep_str = split_list[pep_idx]
            pep_str = pep_str[1:-1]
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
    
def get_pep_from_pepxml_folder(input_folder_str: str, output_file_str: str) -> None:
    fw = open(output_file_str, 'w')
    fw.close()
    pep_set = set()
    onlyfiles = [f for f in listdir(input_folder_str) if isfile(join(input_folder_str, f)) and f.endswith('pep.xml')]
    for input_file_str in onlyfiles:
        input_file_str = join(input_folder_str, input_file_str)
        get_pep_from_pepxml(input_file_str, output_file_str, pep_set)
    
    print('get_pep_from_pepxml_folder is done.')
    
    
def get_pep_from_pepxml(input_file_str: str, output_file_str: str, pep_set: set) -> None:
    fw = open(output_file_str, 'a')
    min_int = 10000
    max_int = 0
    num_protein_int = 0
    with open(input_file_str, 'r') as fr:
        for line_str in fr:
            if "<search_hit" in line_str:
                num_protein_int = 0
                split_list = line_str.strip().split()
                for e in split_list:
                    if "num_tot_proteins" in e:
                        split_2_list = e.split('="')
                        num_protein_int = int(split_2_list[1][:-1])
                        break
                if num_protein_int > 1:
                    for e in split_list:
                        if "peptide=" in e:
                            split_2_list = e.split('="')
                            pep_str = split_2_list[1][:-1]
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
    
def complete_proteins_4_pepxml_folder(input_folder_str: str, pep_sip_file_str: str, output_folder_str: str) -> None:
    onlyfiles = [f for f in listdir(input_folder_str) if isfile(join(input_folder_str, f)) and f.endswith('pep.xml')]
    pep_pro_dict = read_sip(pep_sip_file_str)
    for input_file_str in onlyfiles:
        output_file_str = join(output_folder_str, input_file_str)
        input_file_str = join(input_folder_str, input_file_str)
        complete_proteins_4_pepxml(input_file_str, pep_pro_dict, output_file_str)
    
    print('complete_proteins_4_pepxml_folder is done.')

def read_sip(pep_sip_file_str: str) -> dict:
    pep_pro_dict = {}
    with open(pep_sip_file_str, 'r') as fr:
        for line_str in fr:
            split_list = line_str.strip().split('\t')
            if len(split_list) < 2: # peptide too short or too long not in the sip file
                continue
            pro_list = split_list[1:]
            pro_list = Settings.clean_protein(pro_list)
            pep_str = split_list[0][1:-1]
            if pep_str in pep_pro_dict:
                l = pep_pro_dict[pep_str]
                for pro in pro_list:
                    if pro not in l:
                        l.append(pro)
            else:
                l = []
                pep_pro_dict[pep_str] = l
                for pro in pro_list:
                    if pro not in l:
                        l.append(pro)
    return pep_pro_dict  

def complete_proteins_4_pepxml(input_pepxml_str: str, pep_pro_dict: dict, output_pepxml_str:str) -> None:
    
    fw = open(output_pepxml_str, 'w')
    
    with open(input_pepxml_str, 'r') as fr:
        for line_str in fr:
            fw.write(line_str)
            if "<search_hit" in line_str:
                
                num_protein_int = 0
                split_list = line_str.strip().split()
                for e in split_list:
                    if "num_tot_proteins" in e:
                        split_2_list = e.split('="')
                        num_protein_int = int(split_2_list[1][:-1])
                        break
                if num_protein_int > 1:
                    for e in split_list:
                        if "peptide=" in e:
                            split_2_list = e.split('="')
                            pep_str = split_2_list[1][:-1]
                        if "protein=" in e:
                            split_2_list = e.split('="')
                            pro_str = split_2_list[1][:-1]
                    pro_list = pep_pro_dict[pep_str]
                    for pro in pro_list:
                        if pro != pro_str:
                            fw.write('    <alternative_protein protein="')
                            fw.write(pro)
                            fw.write('"/>\n')
    fw.close()
    
def complete_proteins_and_clean(input_file_str, pep_sip_file_str, output_file_str):
    
    fw = open(output_file_str, 'w')
    pep_pro_dict = {}
    with open(pep_sip_file_str, 'r') as fr:
        for line_str in fr:
            split_list = line_str.strip().split('\t')
            if len(split_list) < 2: # peptide too short or too long not in the sip file
                continue
            pro_list = split_list[1:]
            # pro_list = Settings.clean_protein(pro_list)
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

def complete_proteins_and_clean_pep_psm_files(pep_file_str, psm_file_str, pep_sip_file_str, output_pep_file_str, output_psm_file_str):
    
    fw = open(output_pep_file_str, 'w')
    pep_pro_dict = {}
    with open(pep_sip_file_str, 'r') as fr:
        for line_str in fr:
            split_list = line_str.strip().split('\t')
            if len(split_list) < 2: # peptide too short or too long not in the sip file
                continue
            pro_list = split_list[1:]
            pro_list = Settings.clean_protein(pro_list)
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
    with open(pep_file_str, 'r') as fr:
        header_str = fr.readline()
        fw.write(header_str)
        for line_str in fr:
            split_list = line_str.strip().split('\t')
            pep_str = split_list[2]
            if len(pep_str) < 9 or len(pep_str) > 62:
                continue
            if pep_str in pep_pro_dict:
                pro_list = pep_pro_dict[pep_str]
            else:
                fw.write(line_str)
                continue
            if int(split_list[4]) >= len(pro_list):
                fw.write(line_str)
                continue
            split_list[3] = '{' + ','.join(pro_list) + '}'
            split_list[4] = str(len(pro_list))
            protein_type_int = Settings.get_protein_type(pro_list)
            if protein_type_int <= Settings.LabelFwd:
                split_list[5] = 'T'
            else:
                split_list[6] = 'F'  
            fw.write('\t'.join(split_list))
            fw.write('\n')
    fw.close()
    
    fw = open(output_psm_file_str, 'w')
    with open(psm_file_str, 'r') as fr:
        header_str = fr.readline()
        fw.write(header_str)
        for line_str in fr:
            split_list = line_str.strip().split('\t')
            pep_str = split_list[14]
            if len(pep_str) < 9 or len(pep_str) > 62:
                continue
            if pep_str in pep_pro_dict:
                pro_list = pep_pro_dict[pep_str]
            else:
                fw.write(line_str)
                continue
            if int(split_list[16]) >= len(pro_list):
                fw.write(line_str)
                continue
            split_list[15] = '{' + ','.join(pro_list) + '}'
            split_list[16] = str(len(pro_list))
            protein_type_int = Settings.get_protein_type(pro_list)
            if protein_type_int <= Settings.LabelFwd:
                split_list[17] = 'T'
            else:
                split_list[17] = 'F'  
            fw.write('\t'.join(split_list))
            fw.write('\n')                        
    fw.close
    print('complete_proteins_and_clean_pep_psm_files is done.')
    
def main(args=None):
    
    '''
    # extract peptides from pepxml files
    input_folder_str = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/ecoli_samples/comet/soil/pepxml/"
    output_file_str = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/ecoli_samples/comet/soil/pepxml/pep.txt"
    get_pep_from_pepxml_folder(input_folder_str, output_file_str)
    '''
    
    # complete the proteins for pep xml file
    input_folder_str = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/ecoli_samples/comet/soil/pepxml"
    pep_sip_file_str = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/ecoli_samples/comet/soil/pepxml/pep.sip"
    output_folder_str = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/ecoli_samples/comet/soil/corrected_pepxml"
    complete_proteins_4_pepxml_folder(input_folder_str, pep_sip_file_str, output_folder_str)
    
    '''
    # extract peptides from pin files
    input_file_str = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/ecoli_samples/comet/soil/pin/all_re_index.pin"
    output_file_str = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/ecoli_samples/comet/soil/pin/pep.txt"
    get_pep_from_pin(input_file_str, output_file_str)
    '''
    '''
    # complete the proteins and remove too short or too long peptides, also correct lables
    input_file_str = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/ecoli_samples/comet/soil/pin/all_re_index.pin"
    pep_sip_file_str = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/ecoli_samples/comet/soil/pin/pep.sip"
    output_file_str = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/ecoli_samples/comet/soil/pin/all_corrected.pin"
    complete_proteins_and_clean(input_file_str, pep_sip_file_str, output_file_str)
    '''
    '''
    # extract peptides from pep.txt files
    input_file_str = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/Ecoli/Results/Myrimatch_Comet_MsGF_4_iProphet/soil/fdr_0.02/soil_iprophet_mcm.pep.txt"
    output_file_str = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/Ecoli/Results/Myrimatch_Comet_MsGF_4_iProphet/soil/fdr_0.02/pep_list.sip"
    get_pep_from_pep_txt(input_file_str, output_file_str)
    '''

    '''
    # complete the proteins and remove too short or too long peptides, also correct lables
    pep_file_str = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/Ecoli/Results/Myrimatch_Comet_4_iProphet/soil/fdr_0.02/soil_iprophet.pep.txt"
    psm_file_str = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/Ecoli/Results/Myrimatch_Comet_4_iProphet/soil/fdr_0.02/soil_iprophet.psm.txt"
    pep_sip_file_str = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/Ecoli/Results/Myrimatch_Comet_4_iProphet/soil/fdr_0.02/pep_list.sip"
    output_pep_file_str = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/Ecoli/Results/Myrimatch_Comet_4_iProphet/soil/fdr_0.02/corrected/soil_iprophet.pep.txt"
    output_psm_file_str = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/Ecoli/Results/Myrimatch_Comet_4_iProphet/soil/fdr_0.02/corrected/soil_iprophet.psm.txt"
    complete_proteins_and_clean_pep_psm_files(pep_file_str, psm_file_str, pep_sip_file_str, output_pep_file_str, output_psm_file_str)
    '''
    
if __name__ == '__main__':
    sys.exit(main())
    