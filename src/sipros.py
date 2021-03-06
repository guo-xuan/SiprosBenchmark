'''
Created on Sep 7, 2016

@author: xgo
'''

import getopt, sys, os
import numpy as np
import math
import re
from collections import namedtuple
from sklearn import linear_model

# # Import Sipros package modules
import Settings
import parseconfig

get_file_list_with_ext =  Settings.get_file_list_with_ext

# Some hard-coded parameters
label_train_str = Settings.label_train_str
label_test_str = Settings.label_test_str
label_reserve_str = Settings.label_reserve_str


Test_Fwd_Ratio = 1

mass_window_max_int = 0

#feature_name_list = ['ParentCharge', 'MVH', 'Xcorr', 'WDP', 'ScoreAgreement', 'MassDifferent', 'DeltaRP1', 'DeltaRP2', 'DeltaRP3', 'DeltaRS1', 'DeltaRS2', 'DeltaRS3', 'DiffRP1', 'DiffRP2', 'DiffRP3', 'DiffRS1', 'DiffRS2', 'DiffRS3', 'DiffNorRP1', 'DiffNorRP2', 'DiffNorRP3', 'DiffNorRS1', 'DiffNorRS2', 'DiffNorRS3', 'NMC', 'IPSC', 'OPSC', 'UPSC', 'SPSC', 'pep_psm', 'pro_pep']
#                     0               1      2        3      4                 5                6           7           8           9           10          11         12          13         14         15         16         17         18            19            20            21            22            23            24     25      26      27      28        31          32       33    34         35           36         37              38 
feature_name_list = ['ParentCharge', 'MVH', 'Xcorr', 'WDP', 'ScoreAgreement', 'MassDifferent', 'DeltaRP1', 'DeltaRP2', 'DeltaRP3', 'DeltaRS1', 'DeltaRS2', 'DeltaRS3', 'DiffRP1', 'DiffRP2', 'DiffRP3', 'DiffRS1', 'DiffRS2', 'DiffRS3', 'DiffNorRP1', 'DiffNorRP2', 'DiffNorRP3', 'DiffNorRS1', 'DiffNorRS2', 'DiffNorRS3', 'NMC', 'IPSC', 'OPSC', 'UPSC', 'SPSC', 'MassWindow', 'PPC', 'OPSC_U', 'OPSC_Ma', 'OPSC_D_Mi', 'PSC_U', 'PSC_Ma', 'PSC_D_Mi']
feature_selection_list = [0, 1, 2, 3, 4, 5]
ptm_str = ['~', '!', '@', '>', '<', '%', '^', '&', '*', '(', ')', '/', '$']
ptm_selection_list = [0]
# ptm_selection_list = [0, 2, 3, 4]
# ptm_selection_list = [0, 1, 2, 3, 4, 5, 6, 7, 8]
for x in ptm_selection_list:
    feature_name_list.append(ptm_str[x])

LabelFwd = Settings.LabelFwd
LabelRevTrain = Settings.LabelRevTrain
LabelTest = Settings.LabelTest
LabelRevReserve = Settings.LabelRevReserve

# # Class for PepOutFields object
class PsmFields4(namedtuple('PsmFields',
        ['FileName',  # 0
         'ScanNumber',  # 1
         'ParentCharge',  # 2
         'MeasuredParentMass',  # 3
         'ScanType',  # 4
         'SearchName',  # 5
         'IdentifiedPeptide',  # 6
         'OriginalPeptide',  # 7
         'CalculatedParentMass',  # 8
         'MVH',  # 9
         'Xcorr',  # 10
         'WDP',  # 11
         'ProteinNames',  # 12
         'ScoreAgreement', # 13
         'DeltaRP1', # 14
         'DeltaRP2', 
         'DeltaRP3',
         'DeltaRS1',
         'DeltaRS2',
         'DeltaRS3',
         'DiffRP1',
         'DiffRP2',
         'DiffRP3',
         'DiffRS1',
         'DiffRS2',
         'DiffRS3',
         'DiffNorRP1',
         'DiffNorRP2',
         'DiffNorRP3',
         'DiffNorRS1',
         'DiffNorRS2',
         'DiffNorRS3',
         'RetentionTime',
         'Rank',
         'DeltaP'])): # 33 ,
    def __init__(self):
        self.data = self

class PSM:

    iNumScores = 3
    fNeutronMass = 1.00867108694132 # it is Neutron mass
    pattern = re.compile('[^\w\[\]]')

    def __init__(self, psm_field):
        self.FileName = psm_field.FileName
        self.bFileNameChanged = False
        self.ScanNumber = int(psm_field.ScanNumber)
        self.ParentCharge = int(psm_field.ParentCharge)
        self.ScanType = psm_field.ScanType
        self.SearchName = psm_field.SearchName
        self.lfScores = [float(psm_field.MVH), float(psm_field.Xcorr), float(psm_field.WDP)]
        self.ProteinNames = psm_field.ProteinNames.strip()
        self.ScoreAgreement = int(psm_field.ScoreAgreement)
        self.IdentifiedPeptide = psm_field.IdentifiedPeptide
        s1 = ''.join([char if char.isalnum() else '$' for char in self.IdentifiedPeptide ])
        self.PTMscore = s1.count('$') - 2
        self.OriginalPeptide = psm_field.OriginalPeptide
        self.OriginalPeptide = PSM.pattern.sub('', self.IdentifiedPeptide)
        self.protein_list = []
        self.RealLabel = protein_type(self.ProteinNames, self.protein_list)
        self.lRanks = []
        self.iInnerId = 0
        self.fPredictProbability = 0.0
        self.fMassDiff = 0.0
        self.dM = 0.0
        self.MeasuredParentMass = float(psm_field.MeasuredParentMass)
        self.CalculatedParentMass = float(psm_field.CalculatedParentMass)
        self.iMassWindow = 0
        self.set_mass_diff(self.MeasuredParentMass, self.CalculatedParentMass)
        
        self.score_differential_list = []
        self.sRTime = '-1.000'
        self.fRtMeasured = 0.0
        self.fRtPredict = 0.0
        self.fRtPvalue = 0.0
        self.iLocalRank = 0
        self.DeltaP = 'NA'
        if type(psm_field).__name__ == 'PsmFields3':
            self.score_differential_list.extend(float(i) for i in psm_field[14:-2])
        elif type(psm_field).__name__ == 'PsmFields4':
            self.score_differential_list.extend(float(i) for i in psm_field[14:32])
            self.sRTime = psm_field.RetentionTime
            self.fRtMeasured = float(self.sRTime)
            self.DeltaP = psm_field.DeltaP 
            self.iLocalRank = int(psm_field.Rank)
        else:
            self.score_differential_list.extend(float(i) for i in psm_field[14:])
        
        if len(self.score_differential_list) != 18:
            sys.stderr.write('error score\n')
        
        self.NMC = 0
        self.IPSC = 0
        self.OPSC = 0
        self.UPSC = 0 # unique peptide
        self.SPSC = 0 # shared peptide
        self.NRS = 0
        self.PPC = 0
        self.UPPC = 0
        self.SPPC = 0
        self.feature_list = []
        
        self.ML_feature = []
        self.FDR_feature = []
        self.fdr_product = 0
        
        self.OPSC_UMD = [0, 0, 0]
        
        self.SPSC_UMD = [0, 0, 0]
        
        self.TrainingLabel = 0
        
    def get_feature_list(self):
        del self.feature_list[:]

        self.feature_list.append(self.ParentCharge) # 1: 0
        self.feature_list.extend(self.lfScores) # 2, 3, 4: 1, 2, 3
        self.feature_list.append(self.ScoreAgreement) # 5: 4
        self.feature_list.append(abs(self.fMassDiff)) # 6: 5
        self.feature_list.extend(self.score_differential_list) # 7 - 24: 6 - 23
        self.feature_list.append(self.NMC) # 25: 24
        self.feature_list.append((self.IPSC)) # 26: 25
        # self.OPSC = self.OPSC_UMD[0] + self.OPSC_UMD[1]*0 + self.OPSC_UMD[2]*0 - 1.0
        # self.OPSC = self.OPSC_UMD[0]*0.3333 + self.OPSC_UMD[1]*0.5 + self.OPSC_UMD[2] - 1.0
        self.feature_list.append((self.OPSC)) # 27: 26
        self.feature_list.append((self.UPSC)) # 28: 27
        # self.SPSC = self.SPSC_UMD[0] + self.SPSC_UMD[1]*0 + self.SPSC_UMD[2]*0 - 1.0
        # self.SPSC = self.SPSC_UMD[0]*0.3333 + self.SPSC_UMD[1]*0.5 + self.SPSC_UMD[2] - 1.0
        self.feature_list.append((self.SPSC)) # 29: 28
 
        self.feature_list.append(abs(self.iMassWindow)) # 30: 29
        
        # num replicate spectra
        # self.feature_list.append(self.NRS) # 31: 30
        self.feature_list.append((self.PPC)) # 31: 30
        # self.feature_list.append((self.UPPC)) # 31: 30
        # self.feature_list.append((self.SPPC)) # 32: 31
        # self.feature_list.extend(self.OPSC_UMD) # 32 - 35: 31 - 34
        # self.feature_list.extend(self.SPSC_UMD) # 32 - 35: 31 - 34
        
        for c in ptm_selection_list:
            self.feature_list.append(self.IdentifiedPeptide.count(ptm_str[c])) # 32: 31
        
    
    def set_protein_names(self):
        self.ProteinNames = '{' + ','.join(self.protein_list) + '}'
        
    def add_protein(self, protein_l):
        add_bool = False
        for p in protein_l:
            if p not in self.protein_list:
                add_bool = True
                self.protein_list.append(p)
        
        if add_bool:
            self.set_protein_names()

    def set_feature(self, feature_list):
        del feature_list[:]
        #feature_list.append(self.ParentCharge)
        feature_list.extend(self.lfScores)
        
    def set_mass_diff(self, measured_mass, calculated_mass):
        fDiff = calculated_mass - measured_mass
        fTemp = fDiff
        fCeil = 0
        fDown = 0
        if fDiff >= 0:
            fDiff = fTemp
            fCeil = math.ceil(fTemp)*PSM.fNeutronMass
            fFloor = math.floor(fTemp)*PSM.fNeutronMass
            if fFloor > fTemp:
                fFloor -= PSM.fNeutronMass
            if fCeil - PSM.fNeutronMass > fTemp:
                fCeil -= PSM.fNeutronMass
            if fTemp > fCeil - fTemp:
                fTemp = fCeil - fTemp
            if fDiff > fDiff - fFloor:
                fDiff = abs(fDiff - fFloor)
            if abs(fTemp) < abs(fDiff):
                fDiff = fTemp
                self.dM = -fTemp
            else:
                self.dM = fDiff
        else:
            fCeil = math.ceil(fDiff)*PSM.fNeutronMass
            if fCeil < fDiff:
                fCeil += PSM.fNeutronMass
            fFloor = math.floor(fDiff)*PSM.fNeutronMass
            if fFloor + PSM.fNeutronMass < fDiff:
                fFloor += PSM.fNeutronMass
            fDiff = fTemp
            if abs(fTemp) > fCeil - fTemp:
                fTemp = fCeil - fTemp
            if abs(fDiff) > fDiff - fFloor:
                fDiff = fDiff - fFloor
            fTemp = abs(fTemp)
            fDiff = abs(fDiff)
            if fTemp < fDiff:
                fDiff = fTemp
                self.dM = -fTemp
            else:
                self.dM = fDiff
        self.fMassDiff = fDiff
        
        '''
        MassDiffOriginal = measured_mass - calculated_mass
        MassDiff = MassDiffOriginal
        global mass_window_max_int
        for i in range(-mass_window_max_int, mass_window_max_int):
            if abs(MassDiffOriginal - i*PSM.fNeutronMass) < abs(MassDiff):
                MassDiff = MassDiffOriginal - i*PSM.fNeutronMass
                self.iMassWindow = i
        self.fMassDiff = MassDiff
        '''
        
    def set_fdr_product(self):
        val = 1.0
        for x in self.FDR_feature:
            val *= (1.0 - x)
        val = 1.0 - pow(val, 1.0/float(len(self.FDR_feature)))
        self.fdr_product = val
        
    def clean_protein_name(self):
        self.ProteinNames = ""
        l = []
        for sProtein in self.protein_list:
            sProtein.strip()
            if label_train_str == "":
                if sProtein not in l:
                    l.append(sProtein)
            elif not (sProtein.startswith(label_train_str)):
                if sProtein not in l:
                    l.append(sProtein)
        self.ProteinNames = '{'+','.join(l) + '}'
        self.protein_list = l
        
    def set_real_label(self):
        self.RealLabel = protein_type(self.ProteinNames, self.protein_list)


# # Version control
def get_version():
    return "Sipros Ensemble 1.0.1 (Alpha)"

# # Help message
help_message = '''
Usage:
    python sipros_post_processing.py [options]

Inputs:
    -i PSM.tab
    -c Sipros Ensemble configuration file

Options:
    -h/--help
    -v/--version

Outputs:
    -o output directory
'''

# # Parse options
def parse_options(argv):

    opts, _args = getopt.getopt(argv[1:], "hvVi:c:o:x:",
                                    ["help",
                                     "version",
                                     "input",
                                     "config",
                                     "output"])

    # Default working dir and config file
    input_file = ""
    output_folder = ""
    config_file = ""
    debug_code = ""

    # Basic options
    for option, value in opts:
        if option in ("-h", "--help"):
            print(help_message)
            sys.exit(0)
        if option in ("-v", "-V", "--version"):
            print("sipros_post_processing.py V%s" % (get_version()))
            sys.exit(0)
        if option in ("-i", "--input"):
            input_file = value
        if option in ("-o", "--output"):
            output_folder = value
        if option in ("-c", "--config"):
            config_file = value
        if option in ("-x"):
            debug_code = value
            
    if input_file == "" or output_folder == "":
        print(help_message)
        sys.exit(0)

    output_folder = os.path.join(output_folder, '')

    return (input_file, config_file, output_folder, debug_code)

# # Decoy Reverse Forward protein
def protein_type(protein_sequence, lProtein=None):
    sProteins = protein_sequence.replace('{', '')
    sProteins = sProteins.replace('}', '')
    protein_list = sProteins.split(',')
    if lProtein != None:
        del lProtein[:]
        for sProtein in protein_list:
            sProtein = sProtein.strip()
            if sProtein not in lProtein:
                lProtein.append(sProtein)
                
    
    return Settings.get_protein_type(protein_list)

# # read the psm table
def read_psm_table(input_file):
    
    sip_files_list = []
    
    # check if it is a folder
    if os.path.isdir(input_file) == True:
        lFileList = get_file_list_with_ext(input_file, '.tab')
        for sFileName in lFileList:
            sip_files_list.append(sFileName)
    elif ',' in input_file:
        file_split_list = input_file.split(',')
        for file_str in file_split_list:
            sip_files_list.append(file_str)
    else:
        sip_files_list.append(input_file)
    
    psm_list = []
    
    # read line with csv
    for file_str in sip_files_list:
        
        f = open(file_str, 'r')
        # skip header
        _sHeader = f.readline()
        # get data
        for sLine in f:
            sLine = sLine.strip().split('\t')
            PsmFields_obj = PsmFields4._make(sLine)
            psm_obj = PSM(PsmFields_obj)
            psm_list.append(psm_obj)
        
    i = 0
    for oPsm in psm_list:
        oPsm.iInnerId = i
        i += 1
        
    # sorting all PSM, first on file name, than scan number
    psm_list = sorted(psm_list, key=lambda psm: (psm.FileName, psm.ScanNumber))
        
    return (psm_list, None)

# # Division error handling
divide = Settings.divide
FDR_parameter = 1.0

# # FDR calculator
def FDR_calculator(FP, TP):
    FDR_numerator = float(FP) * float(FDR_parameter)
    FDR_denominator = float(TP)
    FDR_accept = True

    if  FDR_denominator == 0:
        FDR_value = 1.0
        FDR_accept = False
    else:
        FDR_value = divide(FDR_numerator, FDR_denominator)
        FDR_accept = True

    return (FDR_accept, float(FDR_value))

def show_Fdr(psm_list, fdr_float, charge_left_given = -1, charge_right_given = -1):
    
    # list_sorted = sorted(psm_list, key=lambda x: (x.fPredictProbability, 1 - x.fRankProduct) , reverse=True)
    list_sorted = sorted(psm_list, key=lambda x: (-x.fPredictProbability, -x.fMassDiff, -x.PTMscore, x.IdentifiedPeptide))
    Fwd_num = 0
    Rev_num = 0
    Shu_num = 0
    Best_list = [0, 0, 0]

    
    psm_filtered_list = []
    cutoff_probability = 1000.0
    # without considering training label
    for oPsm in list_sorted:
        if charge_left_given != -1 and (oPsm.ParentCharge < charge_left_given or oPsm.ParentCharge > charge_right_given):
            continue
        if oPsm.RealLabel == LabelFwd:
            Fwd_num += 1
        elif oPsm.RealLabel == LabelRevTrain:
            Rev_num += 1
        elif oPsm.RealLabel == LabelTest:
            Shu_num += 1
        else:
            sys.stderr.write('error 768.\n')
        (FDR_accept, FDR_value) = FDR_calculator(Shu_num, Fwd_num)
        if (FDR_accept is True) and (FDR_value <= fdr_float):
            if (Best_list[0] + Best_list[2]) < (Fwd_num + Shu_num):
                Best_list = [Fwd_num, Rev_num, Shu_num]
                cutoff_probability = oPsm.fPredictProbability
            

    for oPsm in list_sorted:
        if charge_left_given != -1 and (oPsm.ParentCharge < charge_left_given or oPsm.ParentCharge > charge_right_given):
            continue
        if oPsm.fPredictProbability >= cutoff_probability:
            psm_filtered_list.append(oPsm)
     
    # sys.stdout.write('\t'+str(len(psm_filtered_list))+'\n')
    # sys.stdout.write("{:,d}\n{:.2f}%\n{:.2f}%\n{:,d}\n{:.2f}%\n".format(Best_list[0], 100.0*float(Best_list[1])/float(Best_list[0]), 100.0*float(Best_list[2])/float(Best_list[0]), best_fwr_pep, 100.0*float(best_shu_pep)/float(best_fwr_pep)))
    # sys.stdout.write("%d\t[%.2f%%]\t[%.2f%%]\t%d\t[%.2f%%]\n" % (Best_list[0], 100.0*float(Best_list[2])/float(Best_list[0]), 100.0*float(Best_list[1])/float(Best_list[0]), best_fwr_pep, 100.0*float(best_shu_pep)/float(best_fwr_pep)))
    '''
    print(str(charge_left_given))
    print(str(cutoff_probability))
    print(str(Best_list[2]))
    print(str(Best_list[0]))
    '''
    return psm_filtered_list

def show_Fdr_Pep(psm_list, fdr_float, charge_left_given = -1, charge_right_given = -1):
    
    list_sorted = sorted(psm_list, key=lambda x: (-x.fPredictProbability, -x.fMassDiff, -x.PTMscore, x.IdentifiedPeptide))
    
    peptide_set = set()
    num_fwr_pep = 0
    num_rev_pep = 0
    num_shu_pep = 0
    best_fwr_pep = 0
    best_shu_pep = 0
    
    psm_filtered_list = []
    cutoff_probability = 1000.0
    # without considering training label
    for oPsm in list_sorted:
        if charge_left_given != -1 and (oPsm.ParentCharge < charge_left_given or oPsm.ParentCharge > charge_right_given):
            continue
        pep_str = oPsm.IdentifiedPeptide + '_' + str(oPsm.ParentCharge)
        if pep_str not in peptide_set:
            if oPsm.RealLabel == LabelFwd:
                num_fwr_pep += 1
                peptide_set.add(pep_str)
            elif oPsm.RealLabel == LabelTest:
                num_shu_pep += 1
                peptide_set.add(pep_str)
            elif oPsm.RealLabel == LabelRevTrain:
                num_rev_pep += 1
                peptide_set.add(pep_str)
            else:
                sys.stderr.write('Error 71341\n')

        (FDR_accept, FDR_value) = FDR_calculator(num_shu_pep, num_fwr_pep)
        if (FDR_accept is True) and (FDR_value <= fdr_float):
            if (num_shu_pep + num_fwr_pep) > (best_shu_pep + best_fwr_pep):
                cutoff_probability = oPsm.fPredictProbability
                best_fwr_pep = num_fwr_pep
                best_shu_pep = num_shu_pep
            

    for oPsm in list_sorted:
        if charge_left_given != -1 and (oPsm.ParentCharge < charge_left_given or oPsm.ParentCharge > charge_right_given):
            continue
        if oPsm.fPredictProbability >= cutoff_probability:
            psm_filtered_list.append(oPsm)
     
    # sys.stdout.write('\t'+str(len(psm_filtered_list))+'\n')
    # sys.stdout.write("%d\t[%.2f%%]\t[%.2f%%]\t%d\t[%.2f%%]\n" % (Best_list[0], 100.0*float(Best_list[2])/float(Best_list[0]), 100.0*float(Best_list[1])/float(Best_list[0]), best_fwr_pep, 100.0*float(best_shu_pep)/float(best_fwr_pep)))
    '''
    print(str(charge_left_given))
    print(str(cutoff_probability))
    print(str(best_shu_pep))
    print(str(best_fwr_pep))
    '''
    return psm_filtered_list


def re_rank(psm_list, consider_charge_bool = False):
    psm_new_list = []
    psm_dict = {}
    if consider_charge_bool :
        for oPsm in psm_list:
            sId = oPsm.FileName + '_' + str(oPsm.ScanNumber) + '_' + str(oPsm.ParentCharge)
            if sId in psm_dict:
                if oPsm.fPredictProbability > psm_dict[sId].fPredictProbability:
                    psm_dict[sId] = oPsm
                elif oPsm.fPredictProbability == psm_dict[sId].fPredictProbability:
                    if abs(oPsm.fMassDiff) < abs(psm_dict[sId].fMassDiff): 
                        psm_dict[sId] = oPsm
                    elif abs(oPsm.fMassDiff) == abs(psm_dict[sId].fMassDiff):
                        # calculate PTM scores
                        if oPsm.PTMscore < psm_dict[sId].PTMscore:
                            psm_dict[sId] = oPsm
                        elif oPsm.PTMscore == psm_dict[sId].PTMscore:
                            if oPsm.IdentifiedPeptide.upper() < psm_dict[sId].IdentifiedPeptide.upper():
                                psm_dict[sId] = oPsm
                            elif oPsm.IdentifiedPeptide.upper() == psm_dict[sId].IdentifiedPeptide.upper():
                                psm_dict[sId].add_protein(oPsm.protein_list)
                    
            else:
                psm_dict[sId] = oPsm
    else :
        for oPsm in psm_list:
            sId = oPsm.FileName + '_' + str(oPsm.ScanNumber)
            if sId in psm_dict:
                if oPsm.fPredictProbability > psm_dict[sId].fPredictProbability:
                    psm_dict[sId] = oPsm
                elif oPsm.fPredictProbability == psm_dict[sId].fPredictProbability:
                    if abs(oPsm.fMassDiff) < abs(psm_dict[sId].fMassDiff): 
                        psm_dict[sId] = oPsm
                    elif abs(oPsm.fMassDiff) == abs(psm_dict[sId].fMassDiff):
                        # calculate PTM scores
                        if oPsm.PTMscore < psm_dict[sId].PTMscore:
                            psm_dict[sId] = oPsm
                        elif oPsm.PTMscore == psm_dict[sId].PTMscore:
                            if oPsm.IdentifiedPeptide.upper() < psm_dict[sId].IdentifiedPeptide.upper():
                                psm_dict[sId] = oPsm
                            elif oPsm.IdentifiedPeptide.upper() == psm_dict[sId].IdentifiedPeptide.upper():
                                psm_dict[sId].add_protein(oPsm.protein_list)
                    
            else:
                psm_dict[sId] = oPsm

    for _key, value in psm_dict.items():
        psm_new_list.append(value)
    '''
    fw = open("/media/xgo/Seagate/Proteomics/Experiments/SIP/ScoreCompare/filtering/SE_50/tmp/psm_temp.txt", 'w')
    
    for oPsm in psm_new_list:
        if oPsm.RealLabel == LabelRevTrain:
            continue 
        oPsm.clean_protein_name()
        fw.write(oPsm.FileName)
        fw.write('\t')
        fw.write(str(oPsm.ScanNumber))
        fw.write('\t')
        fw.write(str(oPsm.ParentCharge))
        fw.write('\t')
        fw.write('%.3f' % oPsm.MeasuredParentMass)
        fw.write('\t')
        fw.write('%.3f' % oPsm.CalculatedParentMass)
        fw.write('\t')
        fw.write('%.3f' % (oPsm.fMassDiff))
        fw.write('\t')
        fw.write('%.3f' % (1000000*(oPsm.fMassDiff)/oPsm.CalculatedParentMass))
        fw.write('\t')
        fw.write(oPsm.ScanType)
        fw.write('\t')
        fw.write(oPsm.SearchName)
        fw.write('\t')
        fw.write('SiprosEnsemble')
        fw.write('\t')
        fw.write(str(oPsm.fPredictProbability))
        fw.write('\t')
        fw.write('NA')
        fw.write('\t')
        fw.write(oPsm.DeltaP)
        fw.write('\t')
        fw.write(oPsm.IdentifiedPeptide)
        fw.write('\t')
        fw.write(oPsm.OriginalPeptide)
        fw.write('\t')
        fw.write(oPsm.ProteinNames)
        fw.write('\t')
        fw.write(str(len(oPsm.protein_list)))
        fw.write('\t')
        if oPsm.RealLabel == LabelFwd:
            fw.write('T')
        else:
            fw.write('F')
        fw.write('\n')
    
    fw.close()
    '''
    return psm_new_list

def cutoff_filtering(psm_list, config_dict=None, fdr_given=None):
    if fdr_given == None:
        fdr_given = float(config_dict[pro_iden_str + FDR_Threshold_str])
    for oPsm in psm_list:
        oPsm.fPredictProbability = oPsm.lfScores[0]
    psm_new_list = re_rank(psm_list, False)
    psm_return_list = []
    charge_set = [[1, 1], [2, 2], [3, 10000]]
    if config_dict[pro_iden_str + FDR_Filtering_str] == 'PSM':
        for e in charge_set:
            psm_filtered_list_local = show_Fdr(psm_new_list, fdr_given, e[0], e[1])
            psm_return_list.extend(psm_filtered_list_local)
    else:
        for e in charge_set:  
            psm_filtered_list_local = show_Fdr_Pep(psm_new_list, fdr_given, e[0], e[1])
            psm_return_list.extend(psm_filtered_list_local)
            
    # psm_return_list = re_rank(psm_return_list, False)
    
    return psm_return_list



def get_num_missed_cleavage_sites(sIdentifiedSeq, sResiduesBeforeCleavage, sResiduesAfterCleavage):
    count_int = 0
    for i in range(len(sIdentifiedSeq) - 1):
        if sIdentifiedSeq[i] in sResiduesBeforeCleavage and sIdentifiedSeq[i + 1] in sResiduesAfterCleavage:
            count_int += 1
    return count_int
                    
def generate_Prophet_features_group(psm_list, config_dict):
    pep_spectra_U_dict = {}
    pep_spectra_M_dict = {}
    pep_spectra_D_dict = {}
    # iLocalRank = 0, 1, 2, 3
    pep_spectra_UMD_list = [pep_spectra_U_dict, pep_spectra_M_dict, pep_spectra_D_dict, pep_spectra_D_dict]
    for oPsm in psm_list:
        if oPsm.OriginalPeptide not in pep_spectra_UMD_list[oPsm.iLocalRank]:
            pep_spectra_UMD_list[oPsm.iLocalRank][oPsm.OriginalPeptide] = 1
        else:
            pep_spectra_UMD_list[oPsm.iLocalRank][oPsm.OriginalPeptide] += 1
    pro_spectra_U_dict = {}
    pro_spectra_M_dict = {}
    pro_spectra_D_dict = {}
    pro_spectra_UMD_list = [pro_spectra_U_dict, pro_spectra_M_dict, pro_spectra_D_dict, pro_spectra_D_dict]
    for oPsm in psm_list:
        pro_list = oPsm.protein_list
        for pro in pro_list:
            if pro not in pro_spectra_UMD_list[oPsm.iLocalRank]:
                pro_spectra_UMD_list[oPsm.iLocalRank][pro] = 1
            else:
                pro_spectra_UMD_list[oPsm.iLocalRank][pro] += 1
    
    for oPsm in psm_list:
        for i in range(3):
            if oPsm.OriginalPeptide in pep_spectra_UMD_list[i]:
                # if oPsm.OPSC_UMD[i] < pep_spectra_UMD_list[i][oPsm.OriginalPeptide]:
                oPsm.OPSC_UMD[i] = pep_spectra_UMD_list[i][oPsm.OriginalPeptide]
            
        pro_list = oPsm.protein_list

        for pro in pro_list:
            for i in range(3):
                if pro in pro_spectra_UMD_list[i]:
                    if oPsm.SPSC_UMD[i] < pro_spectra_UMD_list[i][pro]:
                        oPsm.SPSC_UMD[i] = pro_spectra_UMD_list[i][pro]

def generate_Prophet_features_test(lPsm, config_dict):
    # peptide with PTM dictionary is for IPSC
    peptide_with_modification_dict = {}
    # peptide without PTM dictionary is for OPSC
    peptide_dict = {}
    peptide_protein_dict = {}
    # psm_set = Set()
    for oPsm in lPsm:
        oPsm.NMC = get_num_missed_cleavage_sites(oPsm.OriginalPeptide, 
                                                 config_dict[pep_iden_str + cleave_after_residues_str],
                                                 config_dict[pep_iden_str + cleave_before_residues_str])
        '''
        unique_id_str = oPsm.FileName + '_' + str(oPsm.ScanNumber) + '_' + oPsm.IdentifiedPeptide
        if unique_id_str in psm_set:
            continue
        else:
            psm_set.add(unique_id_str)
        '''            
        if oPsm.IdentifiedPeptide in peptide_with_modification_dict:
            peptide_with_modification_dict[oPsm.IdentifiedPeptide] += 1
        else:
            peptide_with_modification_dict[oPsm.IdentifiedPeptide] = 1
        if oPsm.OriginalPeptide in peptide_protein_dict:
            pro_list = peptide_protein_dict[oPsm.OriginalPeptide]
            for protein in oPsm.protein_list:
                if protein not in pro_list:
                    pro_list.append(protein)
        else:
            pro_list = []
            for pro in oPsm.protein_list:
                if pro not in pro_list:
                    pro_list.append(pro)
            peptide_protein_dict[oPsm.OriginalPeptide] = pro_list
            

    pattern = re.compile('[^\w\[\]]')
    for key, _value in peptide_with_modification_dict.items():
        peptide_str = pattern.sub('', key)
        if peptide_str in peptide_dict:
            # peptide_dict[peptide_str] += 1
            peptide_dict[peptide_str] += _value
        else:
            # peptide_dict[peptide_str] = 1
            peptide_dict[peptide_str] = _value
    
    # # sibling peptides
    pro_unique_dict = {} # number of unique peptide of a protain
    pro_shared_dict = {} # number of shared peptide of a protain
    # debug
    pro_balanced_shared_dict = {}
    # debug
    num_changed = 0
    changed_flag = False
    # psm_set.clear()
    pro_pep_dict = {}
    pro_unique_pep_dict = {}
    pro_shared_pep_dict = {}
    for oPsm in lPsm:
        pro_list = peptide_protein_dict[oPsm.OriginalPeptide]
        '''
        if len(oPsm.protein_list) != len(pro_list):
            print('check 3')
            print(oPsm.protein_list)
            print(pro_list)
        '''    
        changed_flag = False
        print_flag = True
        for protein in pro_list:
            if not protein in oPsm.protein_list:
                changed_flag = True
                if not print_flag:
                    print(pro_list)
                    print(oPsm.protein_list)
                    print_flag = True
                oPsm.protein_list.append(protein)
        if len(oPsm.protein_list) != len(pro_list):
            print('check 4')
            print(oPsm.protein_list)
            print(pro_list)
            exit(1)
            
        if changed_flag:
            oPsm.set_protein_names()
            oPsm.RealLabel = protein_type(oPsm.ProteinNames)
            # print(oPsm.OriginalPeptide)
            num_changed += 1
        '''
        unique_id_str = oPsm.FileName + '_' + str(oPsm.ScanNumber) + '_' + oPsm.IdentifiedPeptide
        if unique_id_str in psm_set:
            continue
        else:
            psm_set.add(unique_id_str)
        '''
        if len(pro_list) > 1:
            num_pro_float = float(len(pro_list))
            for protein in pro_list:
                if protein in pro_shared_dict: 
                    pro_shared_dict[protein] += 1
                else:
                    pro_shared_dict[protein] = 1
                # debug
                if protein in pro_balanced_shared_dict:
                    pro_balanced_shared_dict[protein] += 1.0/num_pro_float
                else:
                    pro_balanced_shared_dict[protein] = 1.0/num_pro_float
                # debug
        else:
            if pro_list[0] in pro_unique_dict:
                pro_unique_dict[pro_list[0]] += 1
            else:
                pro_unique_dict[pro_list[0]] = 1
        
        for pro in pro_list:
            if pro in pro_pep_dict:
                l = pro_pep_dict[pro]
                if oPsm.OriginalPeptide not in l:
                    l.append(oPsm.OriginalPeptide)
            else:
                l = []
                l.append(oPsm.OriginalPeptide)
                pro_pep_dict[pro] = l
        if len(pro_list) == 1:
            if pro_list[0] in pro_unique_pep_dict:
                l = pro_unique_pep_dict[pro_list[0]]
                if oPsm.OriginalPeptide not in l:
                    l.append(oPsm.OriginalPeptide)
            else:
                l = []
                l.append(oPsm.OriginalPeptide)
                pro_unique_pep_dict[pro_list[0]] = l
        else:
            for pro in pro_list:
                if pro in pro_shared_pep_dict:
                    l = pro_shared_pep_dict[pro]
                    if oPsm.OriginalPeptide not in l:
                        l.append(oPsm.OriginalPeptide)
                else:
                    l = []
                    l.append(oPsm.OriginalPeptide)
                    pro_shared_pep_dict[pro] = l
    if num_changed != 0:
        if not print_flag:
            print("num changed %d" % num_changed)
    
    # collect features
    num_unique_per_pro = 0
    num_shared_per_pro = 0
    num_balanced_shared_per_pro = 0
    
    max_unique_per_psm = 0
    max_shared_per_psm = 0
    max_balanced_shared_per_psm = 0
    
    num_unique_per_psm = 0
    num_shared_per_psm = 0
    num_balanced_shared_per_psm = 0
    
    max_per_psm = 0 # linked to protein
    max_balanced_per_psm = 0 # linked to a protein
    
    max_linked_unique_per_psm = 0
    max_linked_shared_per_psm = 0
    max_linked_balanced_unique_per_psm = 0
    max_linked_balanced_shared_per_psm = 0
    
    for oPsm in lPsm:
        oPsm.IPSC = peptide_with_modification_dict[oPsm.IdentifiedPeptide]
        oPsm.OPSC = peptide_dict[oPsm.OriginalPeptide]
        
        max_unique_per_psm = 0
        max_shared_per_psm = 0
        max_balanced_shared_per_psm = 0
    
        num_unique_per_psm = 0
        num_shared_per_psm = 0
        num_balanced_shared_per_psm = 0
        
        max_per_psm = 0 # linked to protein
        max_balanced_per_psm = 0 # linked to a protein
    
        max_linked_unique_per_psm = 0
        max_linked_shared_per_psm = 0
        max_linked_balanced_unique_per_psm = 0
        max_linked_balanced_shared_per_psm = 0
        
        pro_list = peptide_protein_dict[oPsm.OriginalPeptide]
        
        for protein in pro_list:
            if len(pro_pep_dict[protein]) > oPsm.PPC:
                oPsm.PPC = len(pro_pep_dict[protein])
            if len(pro_list) == 1 and len(pro_unique_pep_dict[protein]) > oPsm.UPPC:
                oPsm.UPPC = len(pro_unique_pep_dict[protein])
            if len(pro_list) > 1 and len(pro_shared_pep_dict[protein]) > oPsm.SPPC:
                oPsm.SPPC = len(pro_shared_pep_dict[protein])
            num_unique_per_pro = 0
            num_shared_per_pro = 0
            num_balanced_shared_per_pro = 0
            
            if protein in pro_unique_dict:
                num_unique_per_pro = pro_unique_dict[protein]
                if len(pro_list) == 1:
                    # pass
                    num_unique_per_pro -= 1
                
            if protein in pro_shared_dict:
                num_shared_per_pro = pro_shared_dict[protein]
                if len(pro_list) > 1:
                    # pass
                    num_shared_per_pro -= 1
            
            if protein in pro_balanced_shared_dict:
                num_balanced_shared_per_pro = pro_balanced_shared_dict[protein]
                if len(pro_list) > 1:
                    num_balanced_shared_per_pro -= 1.0/ float(len(pro_list))
            
            if num_unique_per_pro > max_unique_per_psm:
                max_unique_per_psm = num_unique_per_pro
            if num_shared_per_pro > max_shared_per_psm:
                max_shared_per_psm = num_shared_per_pro
            if num_unique_per_pro + num_shared_per_pro > max_per_psm:
                max_per_psm = num_unique_per_pro + num_shared_per_pro
                max_linked_unique_per_psm = num_unique_per_pro
                max_linked_shared_per_psm = num_shared_per_pro
            num_unique_per_psm += num_unique_per_pro
            num_shared_per_psm += num_shared_per_pro
            
            if num_balanced_shared_per_pro > max_balanced_shared_per_psm:
                max_balanced_shared_per_psm = num_balanced_shared_per_pro
            if num_unique_per_pro + num_balanced_shared_per_pro > max_balanced_per_psm:
                max_balanced_per_psm = num_unique_per_pro + num_balanced_shared_per_pro
                max_linked_balanced_unique_per_psm = num_unique_per_pro
                max_linked_balanced_shared_per_psm = num_balanced_shared_per_pro
            num_balanced_shared_per_psm += num_balanced_shared_per_pro
        
        oPsm.UPSC = num_unique_per_psm
        oPsm.SPSC = num_shared_per_psm
        oPsm.SPSC = num_unique_per_psm + num_shared_per_psm
        oPsm.UPSC = max_unique_per_psm
        oPsm.SPSC = max_shared_per_psm
        # oPsm.SPSC = max_shared_per_psm + max_unique_per_psm
        # oPsm.PPC = oPsm.UPPC + oPsm.SPPC
        # oPsm.UPSC = max_linked_unique_per_psm
        # oPsm.SPSC = max_linked_shared_per_psm
        # oPsm.SPSC = max_linked_balanced_shared_per_psm
        if len(oPsm.protein_list) == 1:
            oPsm.UPSC = 1
        else:
            oPsm.UPSC = 0
        oPsm.SPSC = max_linked_unique_per_psm + max_linked_shared_per_psm
        # debug
        # oPsm.SPSC = max_linked_unique_per_psm
        # debug
        # oPsm.SPSC = float(max_linked_unique_per_psm)/float(len(oPsm.protein_list)) + float(max_linked_shared_per_psm)/float(len(oPsm.protein_list))
        # oPsm.SPSC = max_per_psm
        
        # oPsm.UPSC = max_unique_per_psm
        # oPsm.SPSC = max_shared_per_psm
        
# # Exit system with error message
def die(msg=None):
    if msg is not None:
        print >> sys.stderr, msg
        sys.exit(1)

# # Check file exist
def check_file_exist(filename):

    try:
        with open(filename) as _f: pass
    except IOError as _e:
        print >> sys.stderr, '\nCannot open', filename
        die("Program exit!")

# defaul value
decoy_prefix = 'Rev_'
train_decoy_prefix = 'Rev_1_'
test_decoy_prefix = 'Rev_2_'
min_peptide_per_protein = 2
min_unique_peptide_per_protein = 1
remove_decoy_identification = 'No'

pep_iden_str = '[Peptide_Identification]'
search_type_str = 'Search_Type'
fasta_database_str = 'FASTA_Database'
pro_iden_str = '[Protein_Identification]'
FDR_Threshold_str = 'FDR_Threshold'
training_decoy_prefix_str = 'Training_Decoy_Prefix'
testing_decoy_prefix_str = 'Testing_Decoy_Prefix'
sip_iden_str = '[Stable_Isotope_Probing]'
decoy_prefix_str = 'Decoy_Prefix'
reserved_decoy_prefix_str = 'Reserved_Decoy_Prefix'
FDR_Filtering_str = 'FDR_Filtering'
min_peptide_per_protein_str = 'Min_Peptide_Per_Protein'
min_unique_peptide_per_protein_str = 'Min_Unique_Peptide_Per_Protein'
remove_decoy_identification_str = 'Remove_Decoy_Identification'
cleave_after_residues_str = 'Cleave_After_Residues'
cleave_before_residues_str = 'Cleave_Before_Residues'
Mass_Tolerance_Parent_Ion_str = 'Mass_Tolerance_Parent_Ion'
Parent_Mass_Windows_str = 'Parent_Mass_Windows'
Mass_Tolerance_Parent_Ion_str = 'Mass_Tolerance_Parent_Ion'

## Parse config file
def parse_config(config_filename):

    # Save config values to dictionary
    config_dict = {}    # initialize dictionary

    # Call Yinfeng's parseconfig.py module
    check_file_exist(config_filename)
    
    # Save all config values to dictionary
    all_config_dict = parseconfig.parseConfigKeyValues(config_filename)
    
    # get the mass window
    line_str = all_config_dict[pep_iden_str + Parent_Mass_Windows_str]
    e_list = line_str.strip().split(',')
    global mass_window_max_int
    mass_window_max_int = 0
    for e in e_list:
        v = abs(int(e.strip()))
        if v > mass_window_max_int:
            mass_window_max_int = v
    line_str = all_config_dict[pep_iden_str + Mass_Tolerance_Parent_Ion_str]
    if float(line_str.strip()) > 0.05:
        mass_window_max_int += 1
    
    # return config dictionary
    return all_config_dict

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
        
    def add(self, oPsm):
        self.SpectralCount += 1
        if self.BestScore < oPsm.fPredictProbability:
            self.BestScore = oPsm.fPredictProbability
        self.PSMs.append(oPsm.FileName+'['+str(oPsm.ScanNumber) +']')
        self.ScanType.append(oPsm.ScanType)
        self.SearchName.append(oPsm.SearchName)
        if oPsm.RealLabel == LabelFwd:
            self.TargetMatch = 'T'
        
    def set(self, oPsm):
        self.IdentifiedPeptide = oPsm.IdentifiedPeptide
        self.ParentCharge = oPsm.ParentCharge
        self.OriginalPeptide = oPsm.OriginalPeptide
        self.ProteinNames = oPsm.ProteinNames
        self.ProteinCount = len(oPsm.protein_list)
        self.SpectralCount = 1
        self.BestScore = oPsm.fPredictProbability
        self.PSMs.append(oPsm.FileName+'['+str(oPsm.ScanNumber) +']')
        self.ScanType.append(oPsm.ScanType)
        self.SearchName.append(oPsm.SearchName)
        if oPsm.RealLabel == LabelFwd:
            self.TargetMatch = 'T'
        else:
            self.TargetMatch = 'F'
        
    def __repr__(self):
        l = [self.IdentifiedPeptide,
             str(self.ParentCharge),
             self.OriginalPeptide,
             self.ProteinNames,
             str(self.ProteinCount),
             self.TargetMatch,
             str(self.SpectralCount),
             str(self.BestScore),
             ('{'+','.join(self.PSMs)+'}'),
             ('{'+','.join(self.ScanType)+'}'),
             ('{'+','.join(self.SearchName)+'}')]
        
        return '\t'.join(l) 


def generate_psm_pep_txt(base_out, out_folder, psm_filtered_list, config_dict):
    
    # protein_identification message
    pro_iden_msg = ""
    pro_iden_msg += "#\t########################\n"
    pro_iden_msg += "#\t# Filtering Parameters #\n"
    pro_iden_msg += "#\t########################\n"
    pro_iden_msg += "#\t\n"
    pro_iden_msg += "#\t[Protein_Identification]\n"
    pro_iden_msg += "#\t\n"
    pro_iden_msg += "#\t# the prefix of training decoy sequences' locus IDs in the database\n"
    pro_iden_msg += "#\t" + training_decoy_prefix_str + " = " + str(label_train_str) + "\n"
    pro_iden_msg += "#\t\n"
    pro_iden_msg += "#\t# the prefix of test decoy sequences' locus IDs in the database\n"
    pro_iden_msg += "#\t" + testing_decoy_prefix_str + " = " + str(label_test_str) + "\n"
    pro_iden_msg += "#\t\n"
    pro_iden_msg += "#\t# Level of FDR filtering. Options: \"PSM\" and \"Peptide\"\n"
    pro_iden_msg += "#\t" + FDR_Filtering_str + " = " + config_dict[pro_iden_str + FDR_Filtering_str] + "\n"
    pro_iden_msg += "#\t\n"
    pro_iden_msg += "#\t# FDR threshold for filtering peptide identifications\n"
    pro_iden_msg += "#\t" + FDR_Threshold_str + " = " + config_dict[pro_iden_str + FDR_Threshold_str] + "\n"
    pro_iden_msg += "#\t\n"

    # PSM column message
    psm_column_msg = ""
    psm_column_msg += "#\t################\n"
    psm_column_msg += "#\t# Column Names #\n"
    psm_column_msg += "#\t################\n"
    psm_column_msg += "#\t\n"
    psm_column_msg += "#\t[Column Names]\n"
    psm_column_msg += "#\t\n"
    psm_column_msg += "#\tFilename = Filename of input FT2 file\n"
    psm_column_msg += "#\tScanNumber = Scan number of the PSM\n"
    psm_column_msg += "#\tParentCharge = Charge state of the PSM\n"
    psm_column_msg += "#\tMeasuredParentMass = Measured parent mass\n"
    psm_column_msg += "#\tCalculatedParentMass = Calculated parent mass from peptide sequence\n"
    psm_column_msg += "#\tMassErrorDa = Mass error in Da with 1-Da error correction\n"
    psm_column_msg += "#\tMassErrorPPM = Mass error in PPM with 1-Da error correction\n"
    psm_column_msg += "#\tScanType = Scan type of the PSM\n"
    psm_column_msg += "#\tSearchName = Sipros search name\n"
    psm_column_msg += "#\tScoringFunction = Scoring function used in the search\n"
    psm_column_msg += "#\tScore = Predicted Probability of being true PSM\n"
    psm_column_msg += "#\tDeltaZ = Difference between the best PSM score and the next best PSM of this scan\n"
    psm_column_msg += "#\tDeltaP = Difference between the best modified PSM and its PTM isoform\n"
    psm_column_msg += "#\tIdentifiedPeptide = Identified peptide sequence with potential PTMs and mutations\n"
    psm_column_msg += "#\tOriginalPeptide = Original peptide sequence in the FASTA file\n"
    psm_column_msg += "#\tProteinNames = Names of proteins of the peptide\n"
    psm_column_msg += "#\tProteinCount = Number of proteins that the peptide can be assigned to\n"
    psm_column_msg += "#\tTargetMatch = T for target match and F for decoy match\n"
    psm_column_msg += "#\t\n"

    # pep column message
    pep_column_msg = ""
    pep_column_msg += "#\t################\n"
    pep_column_msg += "#\t# Column Names #\n"
    pep_column_msg += "#\t################\n"
    pep_column_msg += "#\t\n"
    pep_column_msg += "#\t[Column Names]\n"
    pep_column_msg += "#\t\n"
    pep_column_msg += "#\tIdentifiedPeptide = Identified peptide sequence with potential PTMs and mutations\n"
    pep_column_msg += "#\tParentCharge = Charge state of identified peptide\n"
    pep_column_msg += "#\tOriginalPeptide = Original peptide sequence in the FASTA file\n"
    pep_column_msg += "#\tProteinNames = Names of proteins of the peptide\n"
    pep_column_msg += "#\tProteinCount = Number of proteins that the peptide can be assigned to\n"
    pep_column_msg += "#\tTargetMatch = T for target match and F for decoy match\n"
    pep_column_msg += "#\tSpectralCount = Number of PSMs in which the peptide is identified\n"
    pep_column_msg += "#\tBestScore = The best score of those PSMs\n"
    pep_column_msg += "#\tPSMs = List of PSMs for the peptide: FT2_Filename[Scan_Number]\n"
    pep_column_msg += "#\tScanType = Scan type of those PSMs\n"
    pep_column_msg += "#\tSearchName = Sipros search name\n"
    pep_column_msg += "#\t\n"

    # get the FDR, # target, # decoy
    psm_target_int = 0
    psm_decoy_int = 0
    pep_target_int = 0
    pep_decoy_int = 0
    pep_set = set()
    for oPsm in psm_filtered_list:
        if oPsm.RealLabel == LabelRevTrain:
            continue
        pep_str = oPsm.IdentifiedPeptide + '_' + str(oPsm.ParentCharge)
        if pep_str not in pep_set:
            if oPsm.RealLabel == LabelFwd:
                pep_target_int += 1
            else:
                pep_decoy_int += 1
            pep_set.add(pep_str)
        
        if oPsm.RealLabel == LabelFwd:
            psm_target_int += 1
        else:
            psm_decoy_int += 1
    
    # write out into files
    base_out = out_folder + '/' + base_out
    psm_txt_file_str = base_out + ".psm.txt"
    with open(psm_txt_file_str, 'w') as fw:
        fw.write(pro_iden_msg)
        fw.write('#\t########################################\n')
        fw.write('#\t####### PSM Filtering by Sipros ########\n')
        fw.write('#\t########################################\n')
        fw.write('#\t\n')
        fw.write('#\t#######################\n')
        fw.write('#\t# Statistical Results #\n')
        fw.write('#\t#######################\n')
        fw.write('#\t\n')
        fw.write('#\t[Statistical_Results]\n')
        fw.write('#\t\n')
        fw.write('#\t# Numbers of psm after filtering\n')
        fw.write('#\tDecoy_PSMs_After_Filtering = %d\n' % psm_decoy_int)
        fw.write('#\tTarget_PSMs_After_Filtering = %d\n' % psm_target_int)
        fw.write('#\t# PSM FDR = Decoy_PSMs_After_Filtering / Target_PSMs_After_Filtering\n')
        fw.write('#\tPSM_FDR = %.2f%%\n' % ((1/Test_Fwd_Ratio) * 100.0 * psm_decoy_int/psm_target_int))
        fw.write('#\t\n')
        fw.write(psm_column_msg)
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
        for oPsm in psm_filtered_list:
            if oPsm.RealLabel == LabelRevTrain:
                continue 
            oPsm.clean_protein_name()
            fw.write(oPsm.FileName)
            fw.write('\t')
            fw.write(str(oPsm.ScanNumber))
            fw.write('\t')
            fw.write(str(oPsm.ParentCharge))
            fw.write('\t')
            fw.write('%.3f' % oPsm.MeasuredParentMass)
            fw.write('\t')
            fw.write('%.3f' % oPsm.CalculatedParentMass)
            fw.write('\t')
            fw.write('%.3f' % (oPsm.fMassDiff))
            fw.write('\t')
            fw.write('%.3f' % (1000000*(oPsm.fMassDiff)/oPsm.CalculatedParentMass))
            fw.write('\t')
            fw.write(oPsm.ScanType)
            fw.write('\t')
            fw.write(oPsm.SearchName)
            fw.write('\t')
            fw.write('SiprosEnsemble')
            fw.write('\t')
            fw.write(str(oPsm.fPredictProbability))
            fw.write('\t')
            fw.write('NA')
            fw.write('\t')
            fw.write(oPsm.DeltaP)
            fw.write('\t')
            fw.write(oPsm.IdentifiedPeptide)
            fw.write('\t')
            fw.write(oPsm.OriginalPeptide)
            fw.write('\t')
            fw.write(oPsm.ProteinNames)
            fw.write('\t')
            fw.write(str(len(oPsm.protein_list)))
            fw.write('\t')
            if oPsm.RealLabel == LabelFwd:
                fw.write('T')
            else:
                fw.write('F')
            fw.write('\n')
            
    # pep_sub_dict for preparing pep_out
    pep_sub_dict = {}    # initialize dict of list
    for oPsm in psm_filtered_list:
        if oPsm.RealLabel == LabelRevTrain:
            continue
        pep_ID = oPsm.IdentifiedPeptide + '_+_' + str(oPsm.ParentCharge)
        if pep_ID in pep_sub_dict:
            pep_sub_dict[pep_ID].add(oPsm)
        else:
            oPeptide = Peptide()
            oPeptide.set(oPsm)
            pep_sub_dict[pep_ID] = oPeptide
    pep_txt_file_str = base_out + ".pep.txt"
    with open(pep_txt_file_str, 'w') as fw:
        fw.write(pro_iden_msg)
        # statistic results
        fw.write('#\t########################################\n')
        fw.write('#\t####### PSM Filtering by Sipros ########\n')
        fw.write('#\t########################################\n')
        fw.write('#\t\n')
        fw.write('#\t#######################\n')
        fw.write('#\t# Statistical Results #\n')
        fw.write('#\t#######################\n')
        fw.write('#\t\n')
        fw.write('#\t[Statistical_Results]\n')
        fw.write('#\t\n')
        fw.write('#\t# Numbers of peptide after filtering\n')
        fw.write('#\tDecoy_peptides_After_Filtering = %d\n' % pep_decoy_int)
        fw.write('#\tTarget_peptides_After_Filtering = %d\n' % pep_target_int)
        fw.write('#\t# peptide FDR = Decoy_peptides_After_Filtering / Target_peptides_After_Filtering\n')
        fw.write('#\tPeptide_FDR = %.2f%%\n' % ((1/Test_Fwd_Ratio) * 100.0*pep_decoy_int/pep_target_int))
        fw.write('#\t\n')
        fw.write(pep_column_msg)
        # for pep out
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
    
    return (psm_txt_file_str, pep_txt_file_str)

def remark_concensus(psm_list):
    psm_dict = {}
    for oPsm in psm_list:
        unique_id_str = oPsm.FileName + '_' + str(oPsm.ScanNumber) + '_' + oPsm.IdentifiedPeptide
        if unique_id_str in psm_dict:
            psm_dict[unique_id_str] += 1
        else:
            psm_dict[unique_id_str] = 1
    
    psm_set = set()
    for oPsm in psm_list:
        unique_id_str = oPsm.FileName + '_' + str(oPsm.ScanNumber) + '_' + oPsm.IdentifiedPeptide
        count_int = psm_dict[unique_id_str]
        oPsm.iLocalRank = 3 - count_int
        if count_int == 2:
            psm_set.add(oPsm.FileName + '_' + str(oPsm.ScanNumber))
    
    for oPsm in psm_list:
        unique_id_str = oPsm.FileName + '_' + str(oPsm.ScanNumber)
        if unique_id_str in psm_set:
            if oPsm.iLocalRank == 2:
                oPsm.iLocalRank = 3 # Mi

def mass_filter(psm_list, config_dict):
    mass_tolerance = float(config_dict[pro_iden_str + Mass_Tolerance_Parent_Ion_str])
    psm_new_list = []
    for oPsm in psm_list:
        if abs(oPsm.fMassDiff) > mass_tolerance:
            continue
        psm_new_list.append(oPsm)
    return psm_new_list

# find the testing decoy verse forward 
def find_train_test_ratio(config_dict):
    database_str = config_dict['[Peptide_Identification]FASTA_Database']
    num_train_int = 0
    num_test_int = 0
    num_reserved_int = 0
    num_fwd_int = 0
    less_train_str = '>' + label_train_str
    less_test_str = '>' + label_test_str
    less_reserve_str = '>' + label_reserve_str
    with open(database_str, 'r') as f:
        for line_str in f:
            if line_str.startswith('>'):
                if line_str.startswith(less_train_str):
                    num_train_int += 1
                elif line_str.startswith(less_test_str):
                    num_test_int += 1
                elif label_reserve_str!= '' and line_str.startswith(less_reserve_str):
                    num_reserved_int += 1
                else:
                    num_fwd_int += 1
    
    global Test_Fwd_Ratio
    Test_Fwd_Ratio = float(num_test_int) / float(num_fwd_int)
    return Test_Fwd_Ratio

def generatePINPercolator(psm_list, output_str):
    with open(output_str, 'w') as fw:
        fw.write("SpecId\tLabel\tScanNr\tExpMass\tCalcMass\tPEP\tPRO\tMVH\tdeltMVH\tXcorr\tdeltXcorr\tWDP\tDdeltWDP\tPepLen\tdM\tabsdM\tCharge1\tCharge2\tCharge3\tCharge4\tCharge5\tenzInt\tPeptide\tProtein")
        fw.write('\n')
        row_list = []
        for psm in psm_list:
            del row_list[:]
            fileID = psm.FileName[psm.FileName.rfind('_')+1:psm.FileName.rfind('.')]
            row_list.append(str(psm.ScanNumber)+fileID+'_'+str(psm.ParentCharge)+'_1') # SpecId
            if psm.RealLabel == LabelFwd: # Label
                row_list.append("1")
            elif psm.RealLabel == LabelTest:
                row_list.append("1")
            else:
                row_list.append("-1")
            row_list.append(str(psm.ScanNumber)+fileID) # ScanNr
            row_list.append(str(psm.MeasuredParentMass)) # ExpMass
            row_list.append(str(psm.CalculatedParentMass)) # CalcMass
            row_list.append(str(psm.OPSC)) # PEP
            row_list.append(str(psm.SPSC)) # PRO
            row_list.append(str(psm.lfScores[0])) # MVH
            row_list.append(str(psm.score_differential_list[9])) # deltMVH
            row_list.append(str(psm.lfScores[1])) # Xcorr
            row_list.append(str(psm.score_differential_list[10])) # deltXcorr
            row_list.append(str(psm.lfScores[2])) # WDP
            row_list.append(str(psm.score_differential_list[11])) # deltWDP
            row_list.append(str(len(psm.OriginalPeptide) - 2)) # PepLen
            row_list.append(str(psm.dM)) # dM
            if abs(psm.dM) > PSM.fNeutronMass:
                print("error mass difference")
            row_list.append(str(psm.fMassDiff)) # absdM
            for i in range(1, 5): # Charge 1-4
                if psm.ParentCharge == i:
                    row_list.append("1")
                else:
                    row_list.append("0")
            if psm.ParentCharge > 5: # Charge5
                row_list.append("1")
            else:
                row_list.append("0")
            row_list.append(str(psm.NMC)) # enzInt
            row_list.append(psm.IdentifiedPeptide.replace("[", "-.").replace("]", ".-")) # Peptide
            row_list.append("\t".join(psm.protein_list)) # Protein
            fw.write("\t".join(row_list))
            fw.write("\n")

def generatePINPercolator_WDP(psm_list, output_str):
    with open(output_str, 'w') as fw:
        fw.write("SpecId\tLabel\tScanNr\tExpMass\tCalcMass\tWDP\tDdeltWDP\tPEP\tPRO\tPepLen\tdM\tabsdM\tCharge1\tCharge2\tCharge3\tCharge4\tCharge5\tenzInt\tPeptide\tProtein")
        fw.write('\n')
        row_list = []
        for psm in psm_list:
            del row_list[:]
            fileID = psm.FileName[psm.FileName.rfind('_')+1:psm.FileName.rfind('.')]
            row_list.append(str(psm.ScanNumber)+fileID+'_'+str(psm.ParentCharge)+'_1') # SpecId
            if psm.RealLabel == LabelFwd: # Label
                row_list.append("1")
            elif psm.RealLabel == LabelTest:
                row_list.append("1")
            else:
                row_list.append("-1")
            row_list.append(str(psm.ScanNumber)+fileID) # ScanNr
            row_list.append(str(psm.MeasuredParentMass)) # ExpMass
            row_list.append(str(psm.CalculatedParentMass)) # CalcMass
            row_list.append(str(psm.lfScores[2])) # WDP
            row_list.append(str(psm.score_differential_list[11])) # deltWDP
            row_list.append(str(psm.OPSC)) # PEP
            row_list.append(str(psm.SPSC)) # PRO
            row_list.append(str(len(psm.OriginalPeptide) - 2)) # PepLen
            row_list.append(str(psm.dM)) # dM
            if abs(psm.dM) > PSM.fNeutronMass:
                print("error mass difference")
            row_list.append(str(psm.fMassDiff)) # absdM
            for i in range(1, 5): # Charge 1-4
                if psm.ParentCharge > i:
                    row_list.append("0")
                else:
                    row_list.append("1")
            if psm.ParentCharge > 5: # Charge5
                row_list.append("1")
            else:
                row_list.append("0")
            row_list.append(str(psm.NMC)) # enzInt
            row_list.append(psm.IdentifiedPeptide.replace("[", "-.").replace("]", ".-")) # Peptide
            row_list.append("\t".join(psm.protein_list)) # Protein
            fw.write("\t".join(row_list))
            fw.write("\n")


def get_pin_sipros_ensemble(input_tab_file_str, config_file_str, output_file_str):
    
    # get the configuration parameters
    config_dict = parse_config(config_file_str)
    
    # read the big psm table
    (psm_list, _base_out) = read_psm_table(input_tab_file_str)
    
    # get score agreement info
    remark_concensus(psm_list)
    # mass filtering
    psm_list = mass_filter(psm_list, config_dict)
    # generate features
    generate_Prophet_features_test(psm_list, config_dict)
    generate_Prophet_features_group(psm_list, config_dict)
    # set feature all PSMs
    for oPsm in psm_list:
        oPsm.get_feature_list()
    
    generatePINPercolator(psm_list, output_file_str)
    
    print('get_pin_sipros_ensemble is done.')
    
def get_pin_sipros3_wdp_ensemble(input_tab_file_str, config_file_str, output_file_str):
    
    # get the configuration parameters
    config_dict = parse_config(config_file_str)
    
    # read the big psm table
    (psm_list, _base_out) = read_psm_table(input_tab_file_str)
    
    psm_dict = {}
    filename_dict = {}
    idx = 0
    for one_psm in psm_list:
        if one_psm.FileName in filename_dict:
            psm_id_str = filename_dict[one_psm.FileName] + '_' + str(one_psm.ScanNumber)
        else:
            filename_dict[one_psm.FileName] = str(idx)
            psm_id_str = str(idx) + '_' + str(one_psm.ScanNumber)
            idx += 1
        if not psm_id_str in psm_dict:
            psm_dict[psm_id_str] = one_psm
        else:
            if psm_dict[psm_id_str].lfScores[2] < one_psm.lfScores[2]:
                psm_dict[psm_id_str] = one_psm
    
    psm_list = []
    
    for _k, v in psm_dict.items():
        psm_list.append(v)
    
    # get score agreement info
    remark_concensus(psm_list)
    # mass filtering
    psm_list = mass_filter(psm_list, config_dict)
    # generate features
    generate_Prophet_features_test(psm_list, config_dict)
    generate_Prophet_features_group(psm_list, config_dict)
    # set feature all PSMs
    for oPsm in psm_list:
        oPsm.get_feature_list()
    
    generatePINPercolator_WDP(psm_list, output_file_str)
    
    print('get_pin_sipros3_wdp_ensemble is done.')
    
def get_pep_from_psm_table(input_file_str, output_file_str):
    
    pep_set = set()
    
    with open(input_file_str, 'r') as fr:
        header_str = fr.readline()
        header_split_list = header_str.split('\t')
        original_pep_idx = header_split_list.index('OriginalPeptide')
        with open(output_file_str, 'w') as fw:
            for line_str in fr:
                split_list = line_str.strip().split('\t')
                pep_set.add(split_list[original_pep_idx][1:-1])
            for pep in pep_set:
                fw.write(pep)
                fw.write('\n')
    print("get_pep_from_psm_table is done.")
    
def complete_proteins_and_clean(input_file_str, pep_sip_file_str, output_file_str):
    
    pep_pro_dict = {}
    with open(pep_sip_file_str, 'r') as fr:
        for line_str in fr:
            split_list = line_str.strip().split('\t')
            if len(split_list) < 2: # peptide too short or too long not in the sip file
                continue
            pro_list = split_list[1:]
            pro_list = Settings.remove_duplicate_protein(pro_list)
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
        header_split_list = header_str.split('\t')
        original_pep_idx = header_split_list.index('OriginalPeptide')
        protein_idx = header_split_list.index('ProteinNames')
        with open(output_file_str, 'w') as fw:
            fw.write(header_str)
            for line_str in fr:
                split_list = line_str.split('\t')
                pep_str = split_list[original_pep_idx]
                if len(pep_str) < 9 or len(pep_str) > 62:
                    continue
                pro_list = []
                if pep_str in pep_pro_dict:
                    pro_list = pep_pro_dict[pep_str]
                    
                protein_str = split_list[protein_idx]
                protein_original_list = protein_str[1:-1].split(',')
                for one_protein in protein_original_list:
                    if one_protein not in pro_list:
                        pro_list.append(one_protein)
                split_list[protein_idx] = '{' + ','.join(pro_list) + '}'
                
                fw.write('\t'.join(split_list))
    
    print("complete_proteins_and_clean is done.")
         

def main(argv=None):
    '''
    # extract peptides from pin files
    input_file_str = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/Ecoli/Results/Myrimatch_Comet_MsGF_4_SiprosEnsemble/marine/marine.tab"
    output_file_str = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/Ecoli/Results/Myrimatch_Comet_MsGF_4_SiprosEnsemble/marine/marine_pep"
    get_pep_from_psm_table(input_file_str, output_file_str)
    '''
    
    # complete the proteins and remove too short or too long peptides, also correct lables
    input_file_str = "/media/xgo/Seagate/Proteomics/Experiments/BenchmarkRev/Msgf_Comet_Myrimatch_4_Siprosensemble/D10/D10.tab"
    pep_sip_file_str = "/media/xgo/Seagate/Proteomics/Experiments/BenchmarkRev/Msgf_Comet_Myrimatch_4_Siprosensemble/marine.sip"
    output_file_str = "/media/xgo/Seagate/Proteomics/Experiments/BenchmarkRev/Msgf_Comet_Myrimatch_4_Siprosensemble/D10/D10_corrected.tab"
    complete_proteins_and_clean(input_file_str, pep_sip_file_str, output_file_str)
    
    
    # generate sipros ensemble for percolator    
    input_tab_file_str = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/ecoli_samples/sipros_ensemble/soil/EColi_Try_HCD_DE10ppm_CS_1000_NCE30_180_Run.tab"
    config_file_str = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/ecoli_samples/config_files/sipros_ensemble/soil.cfg"
    output_file_str = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/ecoli_samples/sipros_ensemble_4_percolator/soil/all.pin"
    get_pin_sipros_ensemble(input_tab_file_str, config_file_str, output_file_str)
    
    
    '''
    # generate sipros 3 wdp for percolator
    input_tab_file_str = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/ecoli_samples/sipros_ensemble/soil/EColi_Try_HCD_DE10ppm_CS_1000_NCE30_180_Run.tab"
    config_file_str = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/ecoli_samples/config_files/sipros_ensemble/soil.cfg"
    output_file_str = "/media/xgo/Seagate/Proteomics/Experiments/SiprosEnsemble/ecoli_samples/sipros3/soil/all.pin"
    get_pin_sipros3_wdp_ensemble(input_tab_file_str, config_file_str, output_file_str)
    '''
    
if __name__ == '__main__':
    sys.exit(main())
