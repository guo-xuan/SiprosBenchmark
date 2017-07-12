'''
Created on Aug 10, 2016

@author: xgo
'''

import sys, re, os
import Settings

label_train_str = Settings.label_train_str
label_test_str = Settings.label_test_str
label_reserve_str = Settings.label_reserve_str
label_ecoli_str = Settings.label_ecoli_str

LabelEcoli = Settings.LabelEcoli
LabelFwd = Settings.LabelFwd
LabelRevTrain = Settings.LabelRevTrain
LabelRevReserve = Settings.LabelRevReserve
LabelTest = Settings.LabelTest

psm_num_list = []
psm_ecoli_num_list = []
pep_num_list = []
pep_ecoli_num_list = []

# # FDR calculator
FDR_calculator = Settings.FDR_calculator

clean_protein = Settings.clean_protein
        
def get_peptide(pep_str):
    pos_start = pep_str.find('.') + 1
    pos_end = pep_str.rfind('.')
    return pep_str[pos_start:pos_end]

def read_pepXML(input_file, mix_version=False):
    
    psm_dict = {}
    C_pattern = re.compile('C\[160\]')
    M_pattern = re.compile('M\[147\]')
    clean_pattern = re.compile('[">/]')
    scan_id = 0
    charge_id = ''
    original_pep = ''
    identified_pep = ''
    protein_l = []
    iProbability = 0.0
    ntt = 0
    nmc = 0
    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith("<spectrum_query "):
                split_l = line.split(' ')
                for one in split_l:
                    if one.startswith('start_scan='):
                        split_l_2 = one.split('=')
                        scan_id = int(clean_pattern.sub('', split_l_2[-1]))
                    if one.startswith('assumed_charge='):
                        split_l_2 = one.split('=')
                        charge_id = clean_pattern.sub('', split_l_2[-1])
                protein_l = []
                ntt = 2
                nmc = 0
            if line.startswith("<parameter name=\"ntt\""):
                split_l = line.split(' ')
                for one in split_l:
                    if one.startswith('value='):
                        split_l_2 = one.split('=')
                        ntt = int(clean_pattern.sub('', split_l_2[-1]))
            if line.startswith("<parameter name=\"nmc\""):
                split_l = line.split(' ')
                for one in split_l:
                    if one.startswith('value='):
                        split_l_2 = one.split('=')
                        nmc = int(clean_pattern.sub('', split_l_2[-1]))
            if line.startswith("<search_hit"):
                split_l = line.split(' ')
                for one in split_l:
                    if one.startswith('peptide='):
                        split_l_2 = one.split('=')
                        original_pep = clean_pattern.sub('', split_l_2[-1])
                        identified_pep = original_pep
                    if one.startswith("protein="):
                        split_l_2 = one.split('=')
                        protein_l.append(clean_pattern.sub('', split_l_2[-1]))
            if line.startswith("<modification_info modified_peptide"):
                split_l = line.split(' ')
                for one in split_l:
                    if one.startswith('modified_peptide='):
                        split_l_2 = one.split('=')
                        identified_pep = C_pattern.sub('C', (clean_pattern.sub('', split_l_2[-1])))
                        identified_pep = M_pattern.sub('M~', (clean_pattern.sub('', split_l_2[-1])))
            if line.startswith("<alternative_protein"):
                split_l = line.split(' ')
                for one in split_l:
                    if one.startswith('protein='):
                        split_l_2 = one.split('=')
                        tmp_str = clean_pattern.sub('', split_l_2[-1])
                        if tmp_str not in protein_l:
                            protein_l.append(tmp_str)
            if line.startswith("<interprophet_result "):
                split_l = line.split(' ')
                for one in split_l:
                    if one.startswith('probability='):
                        split_l_2 = one.split('=')
                        iProbability = float(clean_pattern.sub('', split_l_2[-1]))            
            if line.startswith("</spectrum_query>"):
                if ntt != 2 or nmc > 3:
                    continue
                if len(original_pep) < 7 or len(original_pep) > 60:
                    continue
                if mix_version:
                    protein_type_int = Settings.get_protein_type_ecoli(protein_l)
                else:
                    protein_type_int = Settings.get_protein_type(protein_l)
                if scan_id in psm_dict:
                    if iProbability > psm_dict[scan_id][1]:
                        psm_dict[scan_id] = (protein_type_int, iProbability, scan_id, charge_id, identified_pep, original_pep, protein_l)
                else:
                    psm_dict[scan_id] = (protein_type_int, iProbability, scan_id, charge_id, identified_pep, original_pep, protein_l)
                protein_l = []
    
    psm_list = []   
    
    for _key, value in psm_dict.items():
        psm_list.append(value)
        
    psm_list_sorted = sorted(psm_list, key = lambda psm: psm[1], reverse = True)
    # psm_list_sorted = sorted(psm_list, key = lambda psm: psm[1])
    # psm_list_sorted = sorted(psm_list, key = lambda psm: psm[3])
    
    return psm_list_sorted

def try_fdr(psm_list_sorted, fdr_f):
    num_test_psm = 0
    num_fwr_psm = 0
    best_fwr_psm = 0
    best_test_psm = 0
    score_cutoff = 0
    peptide_set = set()
    num_fwr_pep = 0
    num_test_pep = 0
    best_fwr_pep = 0
    best_test_pep = 0
    for (protein_type_int, iProbability, _scan_id, charge_id, identified_pep, _original_pep, _protein_l) in psm_list_sorted:
        pep_str = identified_pep + '_' + charge_id
        if pep_str not in peptide_set:
            if protein_type_int == LabelFwd:
                num_fwr_pep += 1
                peptide_set.add(pep_str)
            elif protein_type_int == LabelTest:
                num_test_pep += 1
                peptide_set.add(pep_str)
            
        if protein_type_int == LabelFwd:
            num_fwr_psm += 1
        elif protein_type_int == LabelTest:
            num_test_psm += 1
        (FDR_accept, FDR_value) = FDR_calculator(num_test_psm, num_fwr_psm)
        if (FDR_accept is True) and (FDR_value <= fdr_f) and ((num_fwr_psm + num_test_psm) > (best_fwr_psm + best_test_psm)) :
            best_fwr_psm = num_fwr_psm
            best_test_psm = num_test_psm
            score_cutoff = iProbability

        (FDR_accept, FDR_value) = FDR_calculator(num_test_pep, num_fwr_pep)
        if (FDR_accept is True) and (FDR_value <= fdr_f) and ((num_fwr_pep + num_test_pep) > (best_fwr_pep + best_test_pep)) :
            best_fwr_pep = num_fwr_pep
            best_test_pep = num_test_pep
        
        
    # print "# FWD\t# REV\t# SHU"
    print("{:,d} ({:.2f}%)\t{:,d} ({:.2f}%)".format(best_fwr_psm, (100.0*float(best_test_psm)/float(best_fwr_psm)), best_fwr_pep, (100.0*float(best_test_pep)/float(best_fwr_pep))))
    psm_num_list.append(best_fwr_psm)
    pep_num_list.append(best_fwr_pep)
    print(str(score_cutoff))
    return score_cutoff

def try_fdr_ecoli(psm_list_sorted, fdr_f):
    num_test_psm = 0
    num_fwr_psm = 0
    num_true_fwr_psm = 0
    best_fwr_psm = 0
    best_true_fwr_psm = 0
    best_test_psm = 0
    score_cutoff = 0
    peptide_set = set()
    num_fwr_pep = 0
    num_true_fwr_pep = 0
    num_test_pep = 0
    best_fwr_pep = 0
    best_true_fwr_pep = 0
    best_test_pep = 0
    for (protein_type_int, iProbability, _scan_id, charge_id, identified_pep, _original_pep, _protein_l) in psm_list_sorted:
        pep_str = identified_pep + '_' + charge_id
        if pep_str not in peptide_set:
            if protein_type_int <= LabelFwd:
                num_fwr_pep += 1
                peptide_set.add(pep_str)
                if protein_type_int == LabelEcoli:
                    num_true_fwr_pep += 1
            elif protein_type_int == LabelTest:
                num_test_pep += 1
                peptide_set.add(pep_str)
            
        if protein_type_int <= LabelFwd:
            num_fwr_psm += 1
            if protein_type_int == LabelEcoli:
                num_true_fwr_psm += 1
        elif protein_type_int == LabelTest:
            num_test_psm += 1
        (FDR_accept, FDR_value) = FDR_calculator(num_test_psm, num_fwr_psm)
        if (FDR_accept is True) and (FDR_value <= fdr_f) and ((num_fwr_psm + num_test_psm) > (best_fwr_psm + best_test_psm)) :
            best_fwr_psm = num_fwr_psm
            best_true_fwr_psm = num_true_fwr_psm
            best_test_psm = num_test_psm
            score_cutoff = iProbability

        (FDR_accept, FDR_value) = FDR_calculator(num_test_pep, num_fwr_pep)
        if (FDR_accept is True) and (FDR_value <= fdr_f) and ((num_fwr_pep + num_test_pep) > (best_fwr_pep + best_test_pep)) :
            best_fwr_pep = num_fwr_pep
            best_true_fwr_pep = num_true_fwr_pep
            best_test_pep = num_test_pep
        
        
    # print "# FWD\t# REV\t# SHU"
    print("{:,d} ({:,d}) ({:.2f}%)\t{:,d} ({:,d}) ({:.2f}%)".format(best_fwr_psm, best_true_fwr_psm, (100.0*float(best_test_psm)/float(best_fwr_psm)), best_fwr_pep, best_true_fwr_pep, (100.0*float(best_test_pep)/float(best_fwr_pep))))
    psm_num_list.append(best_fwr_psm)
    pep_num_list.append(best_fwr_pep)
    psm_ecoli_num_list.append(best_true_fwr_psm)
    pep_ecoli_num_list.append(best_true_fwr_pep)
    print(str(score_cutoff))
    return score_cutoff

class Peptide:
    
    def __init__(self):
        self.IdentifiedPeptide = ''
        self.ParentCharge = ''
        self.OriginalPeptide = ''
        self.ProteinNames = []
        self.ProteinCount = 0
        self.SpectralCount = 0
        self.BestScore = 0.0
        self.PSMs = []
        self.ScanType = []
        self.SearchName = []
        
    def add(self, proteinname_list, score_float, filename_str, scannum_str, scantype_str, searchname_str, label_int):
        self.SpectralCount += 1
        if self.BestScore < score_float:
            self.BestScore = score_float
        self.PSMs.append(filename_str +'['+str(scannum_str) +']')
        self.ScanType.append(scantype_str)
        self.SearchName.append(searchname_str)
        if label_int == LabelFwd:
            self.TargetMatch = 'T'
        self.ProteinNames += ',' + ','.join(proteinname_list)
        
        
    def set(self, identified_pep_str, parentcharge_str, original_pep_str, proteinname_list, score_float, filename_str, scannum_str, scantype_str, searchname_str, label_int):
        self.IdentifiedPeptide = identified_pep_str
        self.ParentCharge = parentcharge_str
        self.OriginalPeptide = original_pep_str
        self.ProteinNames = ','.join(proteinname_list)
        self.ProteinCount = 1
        self.SpectralCount = 1
        self.BestScore = score_float
        self.PSMs.append(filename_str +'['+str(scannum_str) +']')
        self.ScanType.append(scantype_str)
        self.SearchName.append(searchname_str)
        if label_int == LabelFwd:
            self.TargetMatch = 'T'
        else:
            self.TargetMatch = 'F'
        
    def __repr__(self):
        protein_list = self.ProteinNames.split(',')
        new_list = []
        for protein_str in protein_list:
            if protein_str not in new_list:
                new_list.append(protein_str)
        l = [self.IdentifiedPeptide,
             str(self.ParentCharge),
             self.OriginalPeptide,
             '{' + ','.join(new_list) + '}',
             str(str(len(new_list))),
             self.TargetMatch,
             str(self.SpectralCount),
             str(self.BestScore),
             ('{'+','.join(self.PSMs)+'}'),
             ('{'+','.join(self.ScanType)+'}'),
             ('{'+','.join(self.SearchName)+'}')]
        
        return '\t'.join(l) 

def generate_psm_pep(input_file_str, psm_list_sorted, score_cutoff, fdr):
    filename_str = input_file_str[input_file_str.rfind('/')+ 1:input_file_str.rfind('.')]
    folder_str = input_file_str[:input_file_str.rfind('/') + 1]
    
    # create a temporary directory
    folder_str = folder_str + 'fdr_' + str(fdr) +'/'
    if not os.path.exists(folder_str):
        os.makedirs(folder_str)
    else:
        for the_file in os.listdir(folder_str):
            file_path = os.path.join(folder_str, the_file)
            try:
                if os.path.isfile(file_path):
                    os.unlink(file_path)
                # elif os.path.isdir(file_path): shutil.rmtree(file_path)
            except Exception as e:
                print(e)
    psm_file_str = folder_str + filename_str + '.psm.txt'
    pep_file_str = folder_str + filename_str + '.pep.txt'
    pep_sub_dict = {}
    with open(psm_file_str, 'w') as fw:
        # for psm out
        psm_out_list = ['Filename',  # 0
                        'ScanNumber',  # 1
                        'ParentCharge',  # 2
                        'MeasuredParentMass',  # 3
                        'CalculatedParentMass',  # 4
                        'MassErrorDa',  # 5 CalculatedParentMass - MeasuredParentMass
                        'MassErrorPPM',  # 6 MassErrorDa / CalculatedParentMass
                        'ScanType',  # 7
                        'SearchName',  # 8
                        'ScoringFunction',  # 9
                        'Score',  # 10
                        'DeltaZ',  # 11 the difference score between the rank 1 and 2
                        'DeltaP',  # 12
                        'IdentifiedPeptide',  # 13
                        'OriginalPeptide',  # 14
                        'ProteinNames',  # 15
                        'ProteinCount',  # 16
                        'TargetMatch']  # 17
        fw.write('\t'.join(psm_out_list) + '\n')
        for (protein_type_int, iProbability, scan_id, charge_id, identified_pep, original_pep, protein_l) in psm_list_sorted:
            if not (protein_type_int <= LabelFwd or protein_type_int == LabelTest):
                continue 
            if iProbability < score_cutoff:
                continue
            if protein_type_int <= LabelFwd :
                TargetMatch_str = 'T'
                label_int = LabelFwd
            else: 
                TargetMatch_str = 'F'
                label_int = LabelTest
            proteinname_list = clean_protein(protein_l)
            filename_str = 'one'
            scannum_str = str(scan_id)
            parentcharge_str = charge_id
            score_str = str(iProbability)
            identified_pep_str = '['+identified_pep+']'
            original_pep_str = '['+original_pep+']'
                
            pep_ID = identified_pep_str + '_' + parentcharge_str
            if pep_ID in pep_sub_dict:
                pep_sub_dict[pep_ID].add(proteinname_list, iProbability, filename_str, scannum_str, 'NA', 'NA', label_int)
            else:
                oPeptide = Peptide()
                oPeptide.set(identified_pep_str, parentcharge_str, original_pep_str, proteinname_list, iProbability, filename_str, scannum_str, 'NA', 'NA', label_int)
                pep_sub_dict[pep_ID] = oPeptide
            fw.write(filename_str) # 0
            fw.write('\t')
            fw.write(scannum_str) # 1
            fw.write('\t')
            fw.write(parentcharge_str) # 2
            fw.write('\t')
            fw.write('NA') # 3
            fw.write('\t')
            fw.write('NA') # 4
            fw.write('\t')
            fw.write('NA') # 5
            fw.write('\t')
            fw.write('NA') # 6
            fw.write('\t')
            fw.write('NA') # 7
            fw.write('\t')
            fw.write('NA') # 8
            fw.write('\t')
            fw.write('NA') # 9
            fw.write('\t')
            fw.write(score_str) # 10
            fw.write('\t')
            fw.write('NA') # 11
            fw.write('\t')
            fw.write('NA') # 12
            fw.write('\t')
            fw.write(identified_pep_str) # 13
            fw.write('\t')
            fw.write(original_pep_str) # 14
            fw.write('\t')
            fw.write('{'+','.join(proteinname_list)+'}') # 15
            fw.write('\t')
            fw.write(str(len(proteinname_list))) # 16
            fw.write('\t')
            fw.write(TargetMatch_str) # 17
            fw.write('\n')
    with open(pep_file_str, 'w') as fw:
        pep_out_list = ['IdentifiedPeptide',    #0
                        'ParentCharge',         #1
                        'OriginalPeptide',      #2
                        'ProteinNames',         #3
                        'ProteinCount',         #4
                        'TargetMatch',          #5
                        'SpectralCount',        #6 number of PSMs matched to this peptide
                        'BestScore',            #7 the highest score of those PSMs
                        'PSMs',                 #8 a list of PSMs matched to this peptide. Use{Filename[ScanNumber],Filename[ScanNumber]} format
                        'ScanType',             #9
                        'SearchName']           #10
        fw.write('\t'.join(pep_out_list) + '\n')
        for _pep_id, oPeptide in pep_sub_dict.items():
            fw.write(repr(oPeptide))
            fw.write('\n')
        
def get_psm_pep_fdr(input_file_str):
    print("# PSM(FDR)\t# PEP(FDR)")
    psm_list_sorted = read_pepXML(input_file_str)
    for fdr in [0.005, 0.01, 0.02]:
        score_cutoff = try_fdr(psm_list_sorted, fdr/3.0)
        generate_psm_pep(input_file_str, psm_list_sorted, score_cutoff, fdr)
        
    for num in psm_num_list:
        sys.stdout.write('{:,d}\t'.format(num))
    sys.stdout.write('\n')
    for num in pep_num_list:
        sys.stdout.write('{:,d}\t'.format(num))
    sys.stdout.write('\n')
    print('get_psm_pep_fdr is Done.')

# ecoli mixed with marine or soil database
def get_psm_pep_fdr_ecoli_mixed(input_file_str):
    print("# PSM(FDR)\t# PEP(FDR)")
    psm_list_sorted = read_pepXML(input_file_str, True)
    for fdr in [0.005, 0.01, 0.02]:
        score_cutoff = try_fdr_ecoli(psm_list_sorted, fdr/3.0)
        generate_psm_pep(input_file_str, psm_list_sorted, score_cutoff, fdr)
        
    for idx, num in enumerate(psm_num_list):
        sys.stdout.write('{:,d} '.format(num))
        sys.stdout.write('({:,d})\t'.format(psm_ecoli_num_list[idx]))
    sys.stdout.write('\n')
    for idx, num in enumerate(pep_num_list):
        sys.stdout.write('{:,d}'.format(num))
        sys.stdout.write('({:,d})\t'.format(pep_ecoli_num_list[idx]))
    sys.stdout.write('\n')
    print('get_psm_pep_fdr_ecoli_mixed is Done.')
    
def main(argv=None):
    
    # get the psm pep fdr using pep xml files
    input_file_str = '/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/Ecoli/Results/Myrimatch_Comet_4_iProphet/ecoli/ecoli_iprophet.pepXML'
    get_psm_pep_fdr(input_file_str)
    
    '''
    input_file_str = '/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/Ecoli/Results/Myrimatch_Comet_4_iProphet/marine/marine_iprophet.pepXML'
    get_psm_pep_fdr_ecoli_mixed(input_file_str)
    '''
    
if __name__ == '__main__':
    sys.exit(main())
