'''
Created on Jun 27, 2017

@author: xgo
'''

import sys, re

import Settings
import comet_myrimatch_msgf

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


def before_iProphet(input_file, out_file):
    
    pep_str = ''
    pos_int = 0
    mass_int = 0
    modified_pep = ''
    pos_mass_dict = {}
    modified_info_str = ''
    with open(out_file, 'w') as fw:
        with open(input_file, 'r') as fr:
            for line_str in fr:
                if line_str.startswith('<search_hit'):
                    split_list = line_str.split(' ')
                    for s in split_list:
                        if s.startswith('peptide='):
                            pep_str = s.split('=')[1][1:-1]
                elif line_str.startswith('<modification_info>'):
                    pos_mass_dict.clear()
                    modified_info_str = ''
                    continue
                elif line_str.startswith('<mod_aminoacid_mass'):
                    modified_info_str += line_str
                    split_list = line_str.split(' ')
                    for s in split_list:
                        if s.startswith('position='):
                            pos_int = int(s.split('=')[1][1:-1])
                        if s.startswith('mass='):
                            mass_int = float(s.split('=')[1][1:-4])
                            pos_mass_dict[pos_int] = '[' + str(int(mass_int)) + ']'
                    continue
                elif line_str.startswith('</modification_info'):
                    if len(pos_mass_dict) > 2:
                        pass
                        # print 'check'
                    fw.write('<modification_info modified_peptide="')
                    s_l = list(pep_str)
                    s_new_l = []
                    for i in range(len(s_l)):
                        s_new_l.append(s_l[i])
                        if (i+1) in pos_mass_dict:
                            if s_l[i] == 'C':
                                continue
                            s_new_l.append(pos_mass_dict[i+1])
                    modified_pep = ''.join(s_new_l) + '">\n'
                    fw.write(modified_pep)
                    fw.write(modified_info_str)
                fw.write(line_str)
                
    print('before_iProphet is done')

import comet_myrimatch_msgf
get_id_pair_from_mzid_4msgf = comet_myrimatch_msgf.get_id_pair_from_mzid_4msgf

# scan number in msgf is messed up, correct it using mzid file
def correct_scan_number(input_xml_str, input_mzid_str, output_xml_str):
    
    id_dict = get_id_pair_from_mzid_4msgf(input_mzid_str)
    fw = open(output_xml_str, 'w')
    with open(input_xml_str, 'r') as fr:
        for line_str in fr:
            if line_str.strip().startswith('<spectrum_query spectrum='):
                begin_int = line_str.find('start_scan="')+12
                end_int = line_str.find('" end_scan')
                old_str = line_str[begin_int:end_int]
                new_str = str(id_dict[old_str])
                line_str = line_str.replace(old_str, new_str)
            fw.write(line_str)
    
    fw.close()
    
    print('correct_scan_number is done: ' + input_xml_str)

def main():
    
    # correct the protein labels
    input_file_str = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/ecoli_samples/msgf/soil_subsample/pin/all.pin"
    output_file_str = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/ecoli_samples/msgf/soil_subsample/pin/all_corrected.pin"
    correct_pin(input_file_str, output_file_str)
    
    
    '''
    # correct the scan number
    n = 20
    input_prefix_str = '/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/ecoli_samples/msgf/soil_subsample/pepxml/EColi_Try_HCD_DE10ppm_CS_1000_NCE30_180_subfile'
    input_suffix_str = '.pepXML'
    mzid_prefix_str = '/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/ecoli_samples/msgf/soil_subsample/mzid/'
    mzid_suffix_str = '_msgfplus.mzid'
    output_prefix_str = '/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/ecoli_samples/msgf/soil_subsample/corrected_pepxml/all_subfile'
    for i in range(n):
        input_file_str = input_prefix_str + str(i) + input_suffix_str
        mzid_file_str = mzid_prefix_str + str(i) + mzid_suffix_str
        output_file_str = output_prefix_str + str(i) + input_suffix_str
        correct_scan_number(input_file_str, mzid_file_str, output_file_str)
    '''
    
    '''
    # clean merged pep xml files
    input_file_str = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/ecoli_samples/msgf/soil_subsample/corrected_pepxml/all.pepXML"
    output_file_str = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/ecoli_samples/msgf/soil_subsample/corrected_pepxml/soil_subsample_msgf.pep.xml"
    clean_merged_pepXML(input_file_str, output_file_str)
    '''
    
    '''
    # before iProphet
    input_file_str = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/ecoli_samples/msgf/soil_subsample/corrected_pepxml/soil_subsample_msgf.pep.xml"
    output_file_str = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/ecoli_samples/msgf/soil_subsample/corrected_pepxml/soil_subsample_msgf_corrected.pep.xml"
    before_iProphet(input_file_str, output_file_str)
    input_file_str = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/ecoli_samples/msgf/marine/corrected_pepxml/marine_msgf.pep.xml"
    output_file_str = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/ecoli_samples/msgf/marine/corrected_pepxml/marine_msgf_corrected.pep.xml"
    # before_iProphet(input_file_str, output_file_str)
    input_file_str = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/ecoli_samples/msgf/soil/corrected_pepxml/soil_msgf.pep.xml"
    output_file_str = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/ecoli_samples/msgf/soil/corrected_pepxml/soil_msgf_corrected.pep.xml"
    # before_iProphet(input_file_str, output_file_str)
    '''
    print('All done.')
    
if __name__ == '__main__':
    sys.exit(main())
