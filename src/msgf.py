'''
Created on Jun 27, 2017

@author: xgo
'''

import sys, re

import Settings

# TestRev should be postive label
# internal missed cleavage site should be less than 4
def correct_pin(input_file_str, output_file_str):
    
    C_empty_re = re.compile('C\[[a-zA-Z]+:[0-9]+\]')
    M_start_re = re.compile('M\[[a-zA-Z]+:[0-9]+\]')
    
    with open(output_file_str, 'w') as fw:
        with open(input_file_str, 'r') as fr:
            header_str = fr.readline()
            fw.write(header_str) # header 1
            header_str = header_str.strip()
            header_list = header_str.split('\t')
            enzInt_idx = header_list.index('enzInt')
            pep_idx = header_list.index('Peptide')
            pro_idx = header_list.index('Proteins')
            fw.write(fr.readline()) # header 2
            for line_str in fr:
                split_list = line_str.split('\t')
                if int(split_list[enzInt_idx]) > 3: # missed cleavage sites
                    continue
                pep_str = C_empty_re.sub('C', split_list[pep_idx])
                pep_str = M_start_re.sub('M', pep_str)
                if len(pep_str) < 11 or len(pep_str) > 64:
                    continue
                protein_list = split_list[pro_idx:]
                protein_type = Settings.get_protein_type(protein_list)
                if protein_type == Settings.LabelRevTrain or protein_type == Settings.LabelRevReserve:
                    split_list[1] = "-1"
                else:
                    split_list[1] = "1"
                fw.write('\t'.join(split_list))
    
    print("Correct Pin Done.")
    
def clean_merged_pepXML(input_file_str, output_file_str):
    hold_text = []
    hold_bool = False
    with open(input_file_str, 'r') as fr:
        with open(output_file_str, 'w') as fw:
            for line_str in fr:
                if "<spectrum_query" in line_str:
                    del hold_text[:]
                    hold_bool = False
                if "</msms_run_summary>" in line_str:
                    hold_bool = True
                if hold_bool:
                    hold_text.append(line_str)
                else:
                    fw.write(line_str)
            
            if len(hold_text) > 0:
                for line_str in hold_text:
                    fw.write(line_str)
    
    print("clean_merged_pepXML is done.")
                

def main():
    '''
    # correct the protein labels
    input_file_str = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/Ecoli/Results/MsGF/soil/all.pin"
    output_file_str = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/Ecoli/Results/MsGF/soil/all_corrected.pin"
    correct_pin(input_file_str, output_file_str)
    '''
    
    
    # clean merged pep xml files
    input_file_str = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/Ecoli/Results/MsGF/ecoli_new/xml/all.pepXML"
    output_file_str = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/Ecoli/Results/MsGF/ecoli_new/xml/ecoli_msgf.pep.xml"
    clean_merged_pepXML(input_file_str, output_file_str)


if __name__ == '__main__':
    sys.exit(main())