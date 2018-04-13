'''
Created on Jun 28, 2017

@author: xgo
'''

import xml.etree.ElementTree as ET
import sys, math

mvh_minimum = 100000
xcorr_minimum = 100000
msgf_minimum = 100000

class peptide:
    
    def __init__(self):
        self.identified_pep = ""
        self.original_pep = ""
        self.calculated_mass = None
        # Negative MS-GF+ E value, logged
        self.msgf_score = None
        self.mvh_score = None
        self.xcorr_score = None
        self.protein_list = []
        self.msgf_diff = 0
        self.mvh_diff = 0
        self.xcorr_diff = 0
        self.local_rank = 0
        
class scan:
    
    def __init__(self):
        self.measured_mass = None
        self.charge = None
        self.scan_number = 0
        self.pep_list = []
        self.msgf_min = 0
        self.mvh_min = 0
        self.xcorr_min = 0
        self.msgf_max = 0
        self.mvh_max = 0
        self.xcorr_max = 0
        self.id_str = ''
    
    def add_pep(self, pep):
        self.pep_list.append(pep)
        
    def update_pep_mvh(self, pep):
        merged_bool = False
        for one_pep in self.pep_list:
            if one_pep.identified_pep == pep.identified_pep:
                one_pep.mvh_score = pep.mvh_score
                merged_bool = True
                for one_protein in pep.protein_list:
                    if one_protein not in one_pep.protein_list:
                        one_pep.protein_list.append(one_protein)
                break
        
        if not merged_bool:
            self.pep_list.append(pep)
        
    def update_pep_xcorr(self, pep):
        merged_bool = False
        for one_pep in self.pep_list:
            if one_pep.identified_pep == pep.identified_pep:
                one_pep.xcorr_score = pep.xcorr_score
                merged_bool = True
                for one_protein in pep.protein_list:
                    if one_protein not in one_pep.protein_list:
                        one_pep.protein_list.append(one_protein)
                break
        
        if not merged_bool:
            self.pep_list.append(pep)
            
    def get_percolator_output(self):
        
        self.find_min()
        self.find_diff()
        
        common_list = ['NA'] # filename
        common_list.append(str(self.scan_number)) # scan number
        common_list.append(self.charge) # parent charge
        common_list.append(self.measured_mass) # measured mass
        common_list.append('NA') # scan type
        common_list.append('NA') # search name
        
        common_list_2 = ['0','0','0','0','0','0','0','0','0']
        common_list_3 = ['0','0','0','0','0','0', 'NA']
        results = []
        
        for pep in self.pep_list:
            if pep.mvh_score >= self.mvh_max and self.mvh_max > mvh_minimum:
                l = []
                l.extend(common_list)
                l.append('[' + pep.identified_pep + ']')
                l.append('[' + pep.original_pep +']')
                l.append(pep.calculated_mass)
                l.append(str(pep.mvh_score))
                l.append(str(pep.xcorr_score))
                l.append(str(pep.msgf_score))
                l.append('{'+','.join(pep.protein_list)+'}')
                l.append('1')
                l.extend(common_list_2)
                l.append(str(pep.mvh_diff))
                l.append(str(pep.xcorr_diff))
                l.append(str(pep.msgf_diff))
                l.extend(common_list_3)
                l.append(str(pep.local_rank))
                l.append('NA')
                results.append(l)
            if pep.xcorr_score >= self.xcorr_max and self.xcorr_max > xcorr_minimum:
                l = []
                l.extend(common_list)
                l.append('[' + pep.identified_pep + ']')
                l.append('[' + pep.original_pep +']')
                l.append(pep.calculated_mass)
                l.append(str(pep.mvh_score))
                l.append(str(pep.xcorr_score))
                l.append(str(pep.msgf_score))
                l.append('{'+','.join(pep.protein_list)+'}')
                l.append('1')
                l.extend(common_list_2)
                l.append(str(pep.mvh_diff))
                l.append(str(pep.xcorr_diff))
                l.append(str(pep.msgf_diff))
                l.extend(common_list_3)
                l.append(str(pep.local_rank))
                l.append('NA')
                results.append(l)
            if pep.msgf_score >= self.msgf_max and self.msgf_max > msgf_minimum:
                l = []
                l.extend(common_list)
                l.append('[' + pep.identified_pep + ']')
                l.append('[' + pep.original_pep +']')
                l.append(pep.calculated_mass)
                l.append(str(pep.mvh_score))
                l.append(str(pep.xcorr_score))
                l.append(str(pep.msgf_score))
                l.append('{'+','.join(pep.protein_list)+'}')
                l.append('1')
                l.extend(common_list_2)
                l.append(str(pep.mvh_diff))
                l.append(str(pep.xcorr_diff))
                l.append(str(pep.msgf_diff))
                l.extend(common_list_3)
                l.append(str(pep.local_rank))
                l.append('NA')
                results.append(l)
        
        return results
    
    def find_diff(self):
        pep_sorted_list = sorted(self.pep_list, key=lambda pep: -pep.mvh_score)
        for i in range(len(pep_sorted_list) - 1):
            pep_sorted_list[i].mvh_diff = pep_sorted_list[i].mvh_score - pep_sorted_list[i+1].mvh_score

        pep_sorted_list = sorted(self.pep_list, key=lambda pep: -pep.xcorr_score)
        for i in range(len(pep_sorted_list) - 1):
            pep_sorted_list[i].xcorr_diff = pep_sorted_list[i].xcorr_score - pep_sorted_list[i+1].xcorr_score
            pep_sorted_list[i].local_rank = i
        if len(pep_sorted_list) == 0:
            print('check')
        
        pep_sorted_list[-1].local_rank = len(pep_sorted_list) - 1
            
        pep_sorted_list = sorted(self.pep_list, key=lambda pep: -pep.msgf_score)
        for i in range(len(pep_sorted_list) - 1):
            pep_sorted_list[i].msgf_diff = pep_sorted_list[i].msgf_score - pep_sorted_list[i+1].msgf_score            
    
    def find_min(self):
        l = []
        for pep in self.pep_list:
            if pep.mvh_score:
                l.append(pep.mvh_score)
        if len(l) < 2:
            self.mvh_min = mvh_minimum
            if len(l) == 1:
                self.mvh_max = l[0]
            else:
                self.mvh_max = mvh_minimum
        else:
            l = sorted(l)
            self.mvh_min = l[0]
            self.mvh_max = l[-1]

        l = []
        for pep in self.pep_list:
            if pep.xcorr_score:
                l.append(pep.xcorr_score)
        if len(l) < 2:
            self.xcorr_min = xcorr_minimum
            if len(l) == 1:
                self.xcorr_max = l[0]
            else:
                self.xcorr_max = xcorr_minimum
        else:
            l = sorted(l)
            self.xcorr_min = l[0]
            self.xcorr_max = l[-1]
            
        l = []
        for pep in self.pep_list:
            if pep.msgf_score:
                l.append(pep.msgf_score)
        if len(l) < 2:
            self.msgf_min = msgf_minimum
            if len(l) == 1:
                self.msgf_max = l[0]
            else:
                self.msgf_max = msgf_minimum
        else:
            l = sorted(l)
            self.msgf_score = l[0]
            self.msgf_max = l[-1]
            
        for pep in self.pep_list:
            if pep.msgf_score is None:
                pep.msgf_score = self.msgf_min
            if pep.mvh_score is None:
                pep.mvh_score = self.mvh_min
            if pep.xcorr_score is None:
                pep.xcorr_score = self.xcorr_min
    
            
def read_one_pep_xml_myrimatch(input_file_str, psm_dict):
    global mvh_minimum
    tree = ET.parse(input_file_str)
    root = tree.getroot()
    msms_run_summary = root.find('{http://regis-web.systemsbiology.net/pepXML}msms_run_summary')
    for spectrum_query in msms_run_summary.iter('{http://regis-web.systemsbiology.net/pepXML}spectrum_query'):
        # scan information
        split_list = spectrum_query.attrib['spectrum'].split('.')
        scan_number = int(split_list[-2]) 
        id_str = str(scan_number) + '_' + spectrum_query.attrib['assumed_charge']
        if id_str in psm_dict:
            one_scan = psm_dict[id_str]
        else:
            one_scan = scan()
            one_scan.measured_mass = spectrum_query.attrib['precursor_neutral_mass']
            one_scan.charge = spectrum_query.attrib['assumed_charge']
            psm_dict[id_str] = one_scan
        for search_result in spectrum_query.iter('{http://regis-web.systemsbiology.net/pepXML}search_result'):
            for search_hit in search_result.iter('{http://regis-web.systemsbiology.net/pepXML}search_hit'):
                pep = peptide()
                # set the original peptide
                pep.original_pep = search_hit.attrib['peptide']
                # get the calculated mass
                pep.calculated_mass = search_hit.attrib['calc_neutral_pep_mass']
                # get the identified peptide
                pos_list = []
                for modification_info in search_hit.findall('{http://regis-web.systemsbiology.net/pepXML}modification_info'):
                    for mod_aminoacid_mass in modification_info.iter('{http://regis-web.systemsbiology.net/pepXML}mod_aminoacid_mass'):
                        pos_int = int(mod_aminoacid_mass.attrib['position'])
                        if pep.original_pep[pos_int-1] == 'M':
                            pos_list.append(pos_int)
                pep.identified_pep = pep.original_pep
                if len(pos_list) > 0:
                    pos_list = sorted(pos_list)
                    for idx, v in enumerate(pos_list):
                        pep.identified_pep = pep.identified_pep[:(v+idx)] + '~' + pep.identified_pep[(v+idx):]
                # get the negative log e valeu score
                for search_score in search_hit.iter('{http://regis-web.systemsbiology.net/pepXML}search_score'):
                    if search_score.attrib['name']=='mvh':
                        pep.mvh_score = float(search_score.attrib['value'])
                        if pep.mvh_score < mvh_minimum:
                            mvh_minimum = pep.mvh_score
                        break
                # get proteins
                pep.protein_list.append(search_hit.attrib['protein'])
                for alternative_protein in search_hit.iter('{http://regis-web.systemsbiology.net/pepXML}alternative_protein'):
                    pep.protein_list.append(alternative_protein.attrib['protein'])
                # add to this scan
                one_scan.update_pep_mvh(pep)
    
    print(input_file_str + ' done.')

def read_one_pep_xml_comet(input_file_str, psm_dict):
    global xcorr_minimum
    tree = ET.parse(input_file_str)
    root = tree.getroot()
    msms_run_summary = root.find('{http://regis-web.systemsbiology.net/pepXML}msms_run_summary')
    for spectrum_query in msms_run_summary.iter('{http://regis-web.systemsbiology.net/pepXML}spectrum_query'):
        # scan information
        split_list = spectrum_query.attrib['spectrum'].split('.')
        scan_number = int(split_list[-2]) 
        id_str = str(scan_number) + '_' + spectrum_query.attrib['assumed_charge']
        if id_str in psm_dict:
            one_scan = psm_dict[id_str]
        else:
            one_scan = scan()
            one_scan.measured_mass = spectrum_query.attrib['precursor_neutral_mass']
            one_scan.charge = spectrum_query.attrib['assumed_charge']
            psm_dict[id_str] = one_scan
        for search_result in spectrum_query.iter('{http://regis-web.systemsbiology.net/pepXML}search_result'):
            for search_hit in search_result.iter('{http://regis-web.systemsbiology.net/pepXML}search_hit'):
                pep = peptide()
                # set the original peptide
                pep.original_pep = search_hit.attrib['peptide']
                # get the calculated mass
                pep.calculated_mass = search_hit.attrib['calc_neutral_pep_mass']
                # get the identified peptide
                pos_list = []
                for modification_info in search_hit.findall('{http://regis-web.systemsbiology.net/pepXML}modification_info'):
                    for mod_aminoacid_mass in modification_info.iter('{http://regis-web.systemsbiology.net/pepXML}mod_aminoacid_mass'):
                        pos_int = int(mod_aminoacid_mass.attrib['position'])
                        if pep.original_pep[pos_int-1] == 'M':
                            pos_list.append(pos_int)
                pep.identified_pep = pep.original_pep
                if len(pos_list) > 0:
                    pos_list = sorted(pos_list)
                    for idx, v in enumerate(pos_list):
                        pep.identified_pep = pep.identified_pep[:(v+idx)] + '~' + pep.identified_pep[(v+idx):]
                # get the negative log e valeu score
                for search_score in search_hit.iter('{http://regis-web.systemsbiology.net/pepXML}search_score'):
                    if search_score.attrib['name']=='xcorr':
                        pep.xcorr_score = float(search_score.attrib['value'])
                        if pep.xcorr_score < xcorr_minimum:
                            xcorr_minimum = pep.xcorr_score
                        break
                # get proteins
                pep.protein_list.append(search_hit.attrib['protein'])
                for alternative_protein in search_hit.iter('{http://regis-web.systemsbiology.net/pepXML}alternative_protein'):
                    pep.protein_list.append(alternative_protein.attrib['protein'])
                # add to this scan
                one_scan.update_pep_xcorr(pep)
    
    print(input_file_str + ' done.')


def read_one_pep_xml_msgf(input_file_str, id_dict, psm_dict):
    global msgf_minimum
    tree = ET.parse(input_file_str)
    root = tree.getroot()
    msms_run_summary = root.find('{http://regis-web.systemsbiology.net/pepXML}msms_run_summary')
    for spectrum_query in msms_run_summary.iter('{http://regis-web.systemsbiology.net/pepXML}spectrum_query'):
        # scan information
        one_scan = scan()
        split_list = spectrum_query.attrib['spectrum'].split('.')
        inner_id_str = split_list[-2]
        one_scan.scan_number = id_dict[inner_id_str]
        one_scan.measured_mass = spectrum_query.attrib['precursor_neutral_mass']
        one_scan.charge = spectrum_query.attrib['assumed_charge']
        one_scan.id_str = str(one_scan.scan_number) + '_' + one_scan.charge
        for search_result in spectrum_query.iter('{http://regis-web.systemsbiology.net/pepXML}search_result'):
            for search_hit in search_result.iter('{http://regis-web.systemsbiology.net/pepXML}search_hit'):
                pep = peptide()
                # set the original peptide
                pep.original_pep = search_hit.attrib['peptide']
                # get the calculated mass
                pep.calculated_mass = search_hit.attrib['calc_neutral_pep_mass']
                # get the identified peptide
                pos_list = []
                for modification_info in search_hit.findall('{http://regis-web.systemsbiology.net/pepXML}modification_info'):
                    for mod_aminoacid_mass in modification_info.iter('{http://regis-web.systemsbiology.net/pepXML}mod_aminoacid_mass'):
                        pos_int = int(mod_aminoacid_mass.attrib['position'])
                        if pep.original_pep[pos_int-1] == 'M':
                            pos_list.append(pos_int)
                pep.identified_pep = pep.original_pep
                if len(pos_list) > 0:
                    pos_list = sorted(pos_list)
                    for idx, v in enumerate(pos_list):
                        pep.identified_pep = pep.identified_pep[:(v+idx)] + '~' + pep.identified_pep[(v+idx):]
                # get the negative log e valeu score
                for search_score in search_hit.iter('{http://regis-web.systemsbiology.net/pepXML}search_score'):
                    if search_score.attrib['name']=='EValue':
                        pep.msgf_score = -math.log(float(search_score.attrib['value']))
                        if pep.msgf_score < msgf_minimum:
                            msgf_minimum = pep.msgf_score
                        break
                # get proteins
                pep.protein_list.append(search_hit.attrib['protein'])
                for alternative_protein in search_hit.iter('{http://regis-web.systemsbiology.net/pepXML}alternative_protein'):
                    pep.protein_list.append(alternative_protein.attrib['protein'])
                # add to this scan
                one_scan.add_pep(pep)
        # add to psm list
        if one_scan.id_str in psm_dict:
            print('check')
        psm_dict[one_scan.id_str] = one_scan
    
    print(input_file_str + ' done.')

def get_id_pair_from_mzid_4msgf(input_file_str):
    id_dict = {}
    inner_id_str = None
    
    with open(input_file_str, 'r') as fr:
        for line_str in fr:
            line_str = line_str.strip()
            if line_str.startswith('<SpectrumIdentificationResult spectrumID='):
                begin_int = line_str.find('index=')+6
                end_int = line_str.find('" spectraData')
                inner_id_str = line_str[begin_int:end_int]
            elif line_str.startswith('<cvParam cvRef="PSI-MS" accession="MS:1001115"'):
                begin_int = line_str.find('value=')+7
                end_int = line_str.find('"/>')
                extra_id_str = int(line_str[begin_int:end_int])
                if inner_id_str:
                    id_dict[inner_id_str] = extra_id_str
                    inner_id_str = None
                else:
                    print('error')
    
    return id_dict

def read_msgf_results(psm_dict):
    
    n = 400
    input_prefix_str = '/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/ecoli_samples/msgf/soil/pepxml/all_subfile'
    input_suffix_str = '.pepXML'
    mzid_prefix_str = '/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/ecoli_samples/msgf/soil/mzid/'
    mzid_suffix_str = '_msgfplus.mzid'
    
    for i in range(n):
        input_file_str = input_prefix_str + str(i) + input_suffix_str
        mzid_file_str = mzid_prefix_str + str(i) + mzid_suffix_str
        id_dict = get_id_pair_from_mzid_4msgf(mzid_file_str)
        read_one_pep_xml_msgf(input_file_str, id_dict, psm_dict)

def read_comet_results(psm_dict):
    n = 20
    input_prefix_str = '/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/ecoli_samples/comet/soil/corrected_pepxml/EColi_Try_HCD_DE10ppm_CS_1000_NCE30_180_subfile'
    input_suffix_str = '.pep.xml'
    
    for i in range(n):
        input_file_str = input_prefix_str + str(i) + input_suffix_str
        read_one_pep_xml_comet(input_file_str, psm_dict)

def read_myrimatch_results(psm_dict):
    n = 20
    input_prefix_str = '/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/ecoli_samples/myrimatch/soil/EColi_Try_HCD_DE10ppm_CS_1000_NCE30_180_subfile'
    input_suffix_str = '.pepXML'
    
    for i in range(n):
        input_file_str = input_prefix_str + str(i) + input_suffix_str
        read_one_pep_xml_myrimatch(input_file_str, psm_dict)

def write_output(psm_list, output_file_str):
    with open(output_file_str, 'w') as fw:
        header_str = 'FileName\tScanNumber\tParentCharge\tMeasuredParentMass\tScanType\tSearchName\tIdentifiedPeptide\tOriginalPeptide\tCalculatedParentMass\tMVH\tXcorr\tWDP\tProteinNames\tScoreAgreement\tDeltaRP1\tDeltaRP2\tDeltaRP3\tDeltaRS1\tDeltaRS2\tDeltaRS3\tDiffRP1\tDiffRP2\tDiffRP3\tDiffRS1\tDiffRS2\tDiffRS3\tDiffNorRP1\tDiffNorRP2\tDiffNorRP3\tDiffNorRS1\tDiffNorRS2\tDiffNorRS3\tRetentionTime\tLocalRank\tDeltaP'
        fw.write(header_str)
        fw.write('\n')
        for _k, v in psm_list.items():
            if v.pep_list:
                ll = v.get_percolator_output()
                for l in ll:
                    fw.write('\t'.join(l))
                    fw.write('\n')
    
    print('write_output is done.')

def main():
    
    psm_dict = {}
    
    # read msgf data first
    read_msgf_results(psm_dict)
    
    # read comet results
    read_comet_results(psm_dict)
    
    # read myrimatch results
    read_myrimatch_results(psm_dict)
    
    # output results
    output_file_str = '/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/ecoli_samples/comet_myrimatch_msgf_4_sipros_ensemble/soil/soil.tab'
    write_output(psm_dict, output_file_str)
    
    print('All done.')

if __name__ == '__main__':
    sys.exit(main())