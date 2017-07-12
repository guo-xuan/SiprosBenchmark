'''
Created on Jun 28, 2017

@author: xgo
'''

import xml.etree.ElementTree as ET
import sys
import Settings
from sipros import protein_type

mvh_minimum = 100000
xcorr_minimum = 100000
msgf_minimum = 100000

class peptide:
    
    def __init__(self):
        self.identified_pep = ""
        self.original_pep = ""
        self.calculated_mass = None
        self.mvh_score = None
        self.mz_score = None
        self.protein_list = []
        self.mvh_diff = 0
        self.mz_diff = 0
        self.local_rank = 0
        
class scan:
    
    def __init__(self):
        self.measured_mass = None
        self.charge = None
        self.scan_number = 0
        self.pep_list = []
        self.mvh_max = 0
    
    def add_pep(self, pep):
        self.pep_list.append(pep)
            
    def get_percolator_output(self):
        
        self.find_diff()
        
        # SpecId\tLabel\tScanNr\tExpMass\tCalcMass\tmvh\tmz\tdeltamvh\tdeltamz\tPepLen\tdM\tabsdM\tCharge1\tCharge2\tCharge3\tCharge4\tCharge5\tenzInt\tPeptide\tProtein
        
        common_list = [str(self.scan_number)+'_'+str(self.charge)+'_1'] # scan number
        common_list.append('1')
        common_list.append(str(self.scan_number)) # scan number
        common_list.append(self.measured_mass) # measured mass
        
        results = []
        
        for pep in self.pep_list:
            if pep.mvh_score >= self.mvh_max:
                l = []
                l.extend(common_list)
                protein_type_int = Settings.get_protein_type(pep.protein_list)
                if protein_type_int == Settings.LabelFwd or protein_type_int == Settings.LabelTest:
                    l[1] = '1'
                else:
                    l[1] = '-1'
                l.append(pep.calculated_mass)
                l.append(str(pep.mvh_score))
                l.append(str(pep.mvh_diff))
                l.append(str(pep.mz_score))
                l.append(str(pep.mz_diff))
                l.append(str(len(pep.original_pep)))
                dM, absdM = Settings.get_mass_diff(float(self.measured_mass), float(pep.calculated_mass))
                l.append(str(dM))
                l.append(str(absdM))
                c = int(self.charge)
                for i in range(1, 5): # Charge 1-4
                    if c == i:
                        l.append("1")
                    else:
                        l.append("0")
                if c > 5: # Charge5
                    l.append("1")
                else:
                    l.append("0")
                l.append(str(Settings.get_num_missed_cleavage_sites(pep.original_pep)))
                l.append('-.' + pep.identified_pep + '.-')
                l.append('\t'.join(pep.protein_list))
                results.append(l)
        
        return results
    
    def find_diff(self):
        pep_sorted_list = sorted(self.pep_list, key=lambda pep: -pep.mvh_score)
        self.mvh_max = pep_sorted_list[0].mvh_score
        for i in range(len(pep_sorted_list) - 1):
            pep_sorted_list[i].mvh_diff = pep_sorted_list[i].mvh_score - pep_sorted_list[i+1].mvh_score

        pep_sorted_list = sorted(self.pep_list, key=lambda pep: -pep.mz_score)
        for i in range(len(pep_sorted_list) - 1):
            pep_sorted_list[i].mz_diff = pep_sorted_list[i].mz_score - pep_sorted_list[i+1].mz_score
      
    
def read_one_pep_xml_myrimatch(input_file_str, psm_list):
    global mvh_minimum
    tree = ET.parse(input_file_str)
    root = tree.getroot()
    msms_run_summary = root.find('{http://regis-web.systemsbiology.net/pepXML}msms_run_summary')
    for spectrum_query in msms_run_summary.iter('{http://regis-web.systemsbiology.net/pepXML}spectrum_query'):
        # scan information
        split_list = spectrum_query.attrib['spectrum'].split('.')
        one_scan = scan()
        one_scan.scan_number = int(split_list[-2])
        one_scan.measured_mass = spectrum_query.attrib['precursor_neutral_mass']
        one_scan.charge = spectrum_query.attrib['assumed_charge']
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
                    if search_score.attrib['name']=='mzFidelity':
                        pep.mz_score = float(search_score.attrib['value'])
                # get proteins
                pep.protein_list.append(search_hit.attrib['protein'])
                for alternative_protein in search_hit.iter('{http://regis-web.systemsbiology.net/pepXML}alternative_protein'):
                    pep.protein_list.append(alternative_protein.attrib['protein'])
                # add to this scan
                one_scan.add_pep(pep)
        if len(one_scan.pep_list) > 0:
            psm_list.append(one_scan)
    
    print(input_file_str + ' done.')


def read_myrimatch_results(psm_dict):
    n = 20
    input_prefix_str = '/media/xgo/Seagate/Proteomics/Experiments/BenchmarkRev/Myrimatch/D10/OSU_D10_FASP_Elite_03202014_subfile'
    input_suffix_str = '.pepXML'
    
    for i in range(n):
        input_file_str = input_prefix_str + str(i) + input_suffix_str
        read_one_pep_xml_myrimatch(input_file_str, psm_dict)

def write_output(psm_list, output_file_str):
    with open(output_file_str, 'w') as fw:
        header_str = 'SpecId\tLabel\tScanNr\tExpMass\tCalcMass\tmvh\tmz\tdeltamvh\tdeltamz\tPepLen\tdM\tabsdM\tCharge1\tCharge2\tCharge3\tCharge4\tCharge5\tenzInt\tPeptide\tProtein'
        fw.write(header_str)
        fw.write('\n')
        for v in psm_list:
            if v.pep_list:
                ll = v.get_percolator_output()
                for l in ll:
                    fw.write('\t'.join(l))
                    fw.write('\n')
    
    print('write_output is done.')

def main():
    
    psm_list = []
    
    # read myrimatch results
    read_myrimatch_results(psm_list)
    
    # output results
    output_file_str = '/media/xgo/Seagate/Proteomics/Experiments/BenchmarkRev/Myrimatch_4_percolator/D10/all.pin'
    write_output(psm_list, output_file_str)
    
    print('All done.')

if __name__ == '__main__':
    sys.exit(main())