'''
Created on Sep 21, 2016

@author: xgo
'''

import sys, re, os
import Settings

FDR_calculator = Settings.FDR_calculator
clean_protein = Settings.clean_protein

psm_num_list = []
pep_num_list = []

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
        if label_int == Settings.LabelFwd:
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
        if label_int == Settings.LabelFwd:
            self.TargetMatch = 'T'
        else:
            self.TargetMatch = 'F'
        
    def __repr__(self):
        protein_list = self.ProteinNames.split(',')
        new_list = []
        for protein_str in protein_list:
            if protein_str not in new_list:
                new_list.append(protein_str)
        self.ProteinCount = len(new_list)
        l = [self.IdentifiedPeptide,
             str(self.ParentCharge),
             self.OriginalPeptide,
             '{' + ','.join(new_list) + '}',
             str(self.ProteinCount),
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
            if not(protein_type_int == Settings.LabelFwd or protein_type_int == Settings.LabelTest):
                continue 
            if iProbability < score_cutoff:
                continue
            if protein_type_int == Settings.LabelFwd :
                TargetMatch_str = 'T'
                label_int = Settings.LabelFwd
            else: 
                TargetMatch_str = 'F'
                label_int = Settings.LabelTest
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


def read_percolator_results(input_file):
    
    psm_list = []
    psm_dict = {}
    identified_pep = ''
    original_pep = ''
    pattern = re.compile('[^a-zA-Z]')
    C_empty_re = re.compile('C\[[a-zA-Z]+:[0-9]+\]')
    M_start_re = re.compile('M\[[a-zA-Z]+:[0-9]+\]')
    with open(input_file, 'r') as f:
        f.readline()
        for line in f:
            word_list = line.strip().split()
            name_list = word_list[0].split('_')
            scan_id = int(name_list[-3])
            charge_id = name_list[-2]
            iProbability = -float(word_list[3])
            identified_pep = word_list[4]
            identified_pep = identified_pep[2:-2]
            identified_pep = M_start_re.sub('M*', identified_pep)
            identified_pep = C_empty_re.sub('C', identified_pep)
            original_pep = pattern.sub('', identified_pep)
            protein_l = word_list[5:]
            protein_type_int = Settings.get_protein_type(protein_l)
            if len(identified_pep) < 7 or len(identified_pep) > 60:
                continue
            if protein_type_int == Settings.LabelRevReserve:
                continue
            if scan_id in psm_dict:
                if iProbability > psm_dict[scan_id][1]:
                    psm_dict[scan_id] = (protein_type_int, iProbability, scan_id, charge_id, identified_pep, original_pep, protein_l)
            else:
                psm_dict[scan_id] = (protein_type_int, iProbability, scan_id, charge_id, identified_pep, original_pep, protein_l)
    
    for _key, value in psm_dict.items():
        psm_list.append(value)
    
    psm_list_sorted = sorted(psm_list, key = lambda psm: psm[1], reverse = True)
    
    return psm_list_sorted

def try_fdr(psm_list_sorted, fdr_f):
    num_rev = 0
    num_fwr = 0
    best_fwr = 0
    best_rev = 0
    score_cutoff = 0
    peptide_set = set()
    num_fwr_pep = 0
    num_rev_pep = 0
    best_fwr_pep = 0
    best_rev_pep = 0
    for (protein_type_int, iProbability, _scan_id, charge_id, identified_pep, _original_pep, _protein_l) in psm_list_sorted:
        pep_charge_str = identified_pep + '_' + charge_id
        if pep_charge_str not in peptide_set:
            if protein_type_int == Settings.LabelFwd:
                num_fwr_pep += 1
                peptide_set.add(pep_charge_str)
            elif protein_type_int == Settings.LabelTest:
                num_rev_pep += 1
                peptide_set.add(pep_charge_str)
            
        if protein_type_int == Settings.LabelFwd:
            num_fwr += 1
        elif protein_type_int == Settings.LabelTest:
            num_rev += 1
        (FDR_accept, FDR_value) = FDR_calculator(num_rev, num_fwr)
        if (FDR_accept is True) and (FDR_value <= fdr_f) and ((num_fwr + num_rev) > (best_fwr + num_rev)) :
            best_fwr = num_fwr
            best_rev = num_rev
            score_cutoff = iProbability
        
        (FDR_accept, FDR_value) = FDR_calculator(num_rev_pep, num_fwr_pep)
        if (FDR_accept is True) and (FDR_value <= fdr_f) and ((num_fwr_pep + num_rev_pep) > (best_fwr_pep + best_rev_pep)) :
            best_fwr_pep = num_fwr_pep
            best_rev_pep = num_rev_pep
        
    # print "# FWD\t# REV\t# SHU"
    print("{:,d} ({:.2f}%)\t{:,d} ({:.2f}%)".format(best_fwr, (100.0*float(best_rev)/float(best_fwr)), best_fwr_pep, (100.0*float(best_rev_pep)/float(best_fwr_pep))))
    psm_num_list.append(best_fwr)
    pep_num_list.append(best_fwr_pep)
    return score_cutoff

    # print("Done.")


def get_fdr_at_3_levels():
    # input_file = '/media/xgo/Seagate/Proteomics/Experiments/Benchmark/Comet/D2/pin/all_complete_protein.txt'
    # input_file = '/media/xgo/Seagate/Proteomics/Experiments/Benchmark/Comet/09272013/pin/all.txt'
    input_file = '/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/Ecoli/Results/MsGF/ecoli_new/all.txt'
    # global shu_str
    # shu_str = 'TestRxx_'
    # input_file = "/media/xgo/Seagate/Proteomics/Experiments/BenchmarkRev/Comet/1002/pin/all_complete_protein.txt"
    # input_file = "/media/xgo/Seagate/Proteomics/Experiments/BenchmarkRev/Comet/D7/pin/all_complete_protein.txt"
    # input_file = '/media/xgo/Seagate/Proteomics/Data/Zhou/Isolife/Results/Comet/pin/run3/all.txt'
    print("# PSM(FDR)\t# PEP(FDR)")
    psm_list_sorted = read_percolator_results(input_file)
    for fdr in [0.005, 0.01, 0.02]:
        score_cutoff = try_fdr(psm_list_sorted, fdr/3.0)
        generate_psm_pep(input_file, psm_list_sorted, score_cutoff, fdr)
    
    for num in psm_num_list:
        sys.stdout.write('{:,d}\t'.format(num))
    sys.stdout.write('\n')
    for num in pep_num_list:
        sys.stdout.write('{:,d}\t'.format(num))
    sys.stdout.write('\n')
    print("get FDR at 3 level is Done.")
    
def main(argv=None):
    
    # get fdr at 0.5 pct, 1 pct, and 2 pct
    get_fdr_at_3_levels()
    
    print("All done.")

if __name__ == '__main__':
    sys.exit(main())
