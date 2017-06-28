'''
Created on Sep 7, 2016

@author: xgo
'''

import getopt, sys, os
import numpy as np
import csv
import math
import re
from sklearn.tree import DecisionTreeClassifier
from collections import namedtuple
from sklearn.ensemble import RandomForestClassifier
from sklearn import linear_model
from sklearn import preprocessing
# from sklearn.neural_network import MLPClassifier
from subprocess import call
from multiprocessing import Process
from multiprocessing import Queue, cpu_count
import scipy.stats as spst

# # Import Sipros package modules
import sipros_post_module
import parseconfig

# # Class for ignoring comments '#' in sipros file
CommentedFile = sipros_post_module.CommentedFile

#feature_name_list = ['ParentCharge', 'MVH', 'Xcorr', 'WDP', 'ScoreAgreement', 'MassDifferent', 'DeltaRP1', 'DeltaRP2', 'DeltaRP3', 'DeltaRS1', 'DeltaRS2', 'DeltaRS3', 'DiffRP1', 'DiffRP2', 'DiffRP3', 'DiffRS1', 'DiffRS2', 'DiffRS3', 'DiffNorRP1', 'DiffNorRP2', 'DiffNorRP3', 'DiffNorRS1', 'DiffNorRS2', 'DiffNorRS3', 'NMC', 'IPSC', 'OPSC', 'UPSC', 'SPSC', 'pep_psm', 'pro_pep']
feature_name_list = ['ParentCharge', 'MVH', 'Xcorr', 'WDP', 'ScoreAgreement', 'MassDifferent', 'DeltaRP1', 'DeltaRP2', 'DeltaRP3', 'DeltaRS1', 'DeltaRS2', 'DeltaRS3', 'DiffRP1', 'DiffRP2', 'DiffRP3', 'DiffRS1', 'DiffRS2', 'DiffRS3', 'DiffNorRP1', 'DiffNorRP2', 'DiffNorRP3', 'DiffNorRS1', 'DiffNorRS2', 'DiffNorRS3', 'NMC', 'IPSC', 'OPSC', 'UPSC', 'SPSC', 'MassWindow', 'PPC', 'OPSC_U', 'OPSC_Ma', 'OPSC_D_Mi', 'PSC_U', 'PSC_Ma', 'PSC_D_Mi']
feature_selection_list = [0, 1, 2, 3, 4, 5]
ptm_str = ['~', '!', '@', '>', '<', '%', '^', '&', '*', '(', ')', '/', '$']
ptm_selection_list = [0]
# ptm_selection_list = [0, 2, 3, 4]
# ptm_selection_list = [0, 1, 2, 3, 4, 5, 6, 7, 8]
for x in ptm_selection_list:
    feature_name_list.append(ptm_str[x])
    

label_train_str = 'Rev_' # 'Rev_'
label_test_str = 'TestRev_' # 'Shu_'
label_reserve_str = 'Rev_2_'
label_ecoli_str = 'lcl'

# label_train_str = 'Rev1_' # 'Rev_'
# label_test_str = 'Rev2_' # 'Shu_'
# label_reserve_str = ''

#label_train_str = 'Sh1_' # 'Rev_'
#label_test_str = 'Sh2_' # 'Shu_'

#label_train_str = 'Rev_'
#label_test_str = 'Dec_'

LabelEcoli = 0
LabelFwd = 1
LabelRevTrain = 2
LabelTest = 3
LabelRevReserve = 4
LabelNotEcoli = 5

mass_tolerance = 0.03

# # Class for PepOutFields object
class PsmFields(namedtuple('PsmFields',
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
         'ScoreAgreement'])):  # 13

    def __init__(self):
        self.data = self
        
# # Class for PepOutFields object
class PsmFields2(namedtuple('PsmFields',
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
         'DeltaRP1',
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
         'DiffNorRS3'])):  

    def __init__(self):
        self.data = self

        
# # Class for PepOutFields object
class PsmFields3(namedtuple('PsmFields',
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
         'DeltaRP1',
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
         'pep_psm',
         'pro_pep'])):  

    def __init__(self):
        self.data = self
'''
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
         'DeltaRP1',
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
         'RetentionTime'])):

    def __init__(self):
        self.data = self
'''

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
         'DeltaRP1',
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

LabelUnknown = 2
LabelPositive = 1
LabelNegative = 0
LabelFiltered = 3

shuffle_num_control = {'3_1': 3, 
                       '3_MX': 31,
                       '3_XW': 11,
                       '3_MW': 10,
                       '3_A': 414,
                       '2_1': 5,
                       '2_MX': 16,
                       '2_XW': 1,
                       '2_MW': 13,
                       '2_A': 671}

'''
shuffle_num_control = {'3_1': 3, 
                       '3_2': 37,
                       '3_3': 414,
                       '2_1': 5,
                       '2_2': 14,
                       '2_3': 671}
'''
feature_list_str = ['ParentCharge', 'MassDiff', 'MVH', 'Xcorr', 'WDP', 'ScoreAgreement', 'DeltaRP1', 'DeltaRP2', 'DeltaRP3', 'DeltaRS1', 'DeltaRS2', 'DeltaRS3', 'DiffRP1', 'DiffRP2', 'DiffRP3', 'DiffRS1', 'DiffRS2', 'DiffRS3', 'DiffNorRP1', 'DiffNorRP2', 'DiffNorRP3', 'DiffNorRS1', 'DiffNorRS2', 'DiffNorRS3']

class PSM:

    iNumScores = 3
    fNeutronMass = 1.00867108694132 # it is Neutron mass
    pattern = re.compile('[^\w\[\]]')

    def __init__(self, psm_field):
        if type(psm_field).__name__ == 'PsmFields3':
            self.pep_psm = int(psm_field.pep_psm)
            self.pro_pep = int(psm_field.pro_pep)
        else:
            self.pep_psm = -1
            self.pro_pep = -1
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
        self.OriginalPeptide = psm_field.OriginalPeptide
        self.OriginalPeptide = PSM.pattern.sub('', self.IdentifiedPeptide)
        self.protein_list = []
        self.RealLabel = protein_type(self.ProteinNames, self.protein_list)
        self.ecoli_label = protein_type_ecoli(self.protein_list)
        self.lRanks = []
        self.TrainingLabel = LabelNegative  # 0: negative 1: positive 2: unknown
        self.fRankProduct = 0.0
        self.iInnerId = 0
        self.fPredictProbability = 0.0
        self.fMassDiff = 0.0
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
            print('error score')
        
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
        
        
    def get_feature_list(self):
        del self.feature_list[:]
        '''
        if self.ParentCharge <= 2:
            self.feature_list.append(2)
        else:
            self.feature_list.append(3)
        '''
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
        '''
        if self.pep_psm != -1:
            self.feature_list.append(self.pep_psm)
            self.feature_list.append(self.pro_pep)
        '''   
        self.feature_list.append(abs(self.iMassWindow)) # 30: 29
        
        # num replicate spectra
        # self.feature_list.append(self.NRS) # 31: 30
        self.feature_list.append((self.PPC)) # 31: 30
        # self.feature_list.append((self.UPPC)) # 31: 30
        # self.feature_list.append((self.SPPC)) # 32: 31
        self.feature_list.extend(self.OPSC_UMD) # 32 - 35: 31 - 34
        self.feature_list.extend(self.SPSC_UMD) # 32 - 35: 31 - 34
        
        for c in ptm_selection_list:
            self.feature_list.append(self.IdentifiedPeptide.count(ptm_str[c])) # 32: 31
        
    
    def set_protein_names(self):
        self.ProteinNames = '{' + ','.join(self.protein_list) + '}'

    def set_feature(self, feature_list):
        del feature_list[:]
        '''
        if self.ParentCharge <= 2:
            feature_list.append(2)
        else:
            feature_list.append(3)
        '''
        #feature_list.append(self.ParentCharge)
        feature_list.extend(self.lfScores)
        '''
        for x in self.lRanks:
            if x == 0:
                feature_list.append(math.log(0.000001))
            else:
                feature_list.append(math.log(x))
        '''
        #feature_list.append(self.lfScores[0])
        #feature_list.append(self.lfScores[1])
        #feature_list.append(self.lfScores[2])
        #feature_list.append(self.ScoreAgreement)
        #feature_list.append(self.fMassDiff)
        
    def set_mass_diff(self, measured_mass, calculated_mass):
        MassDiffOriginal = measured_mass - calculated_mass
        MassDiff = MassDiffOriginal
        for i in range(-4, 4):
            if abs(MassDiffOriginal - i*PSM.fNeutronMass) < abs(MassDiff):
                MassDiff = MassDiffOriginal - i*PSM.fNeutronMass
                self.iMassWindow = i
        self.fMassDiff = MassDiff
        
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
            if not (sProtein.startswith(label_train_str)):
                if sProtein not in l:
                    l.append(sProtein)
        self.ProteinNames = '{'+','.join(l) + '}'
        self.protein_list = l

# # Version control
def get_version():
    return "1.0.1 (Alpha)"

# # Help message
help_message = '''
Usage:
    python xxx.py [options]

Inputs:
    input yyy
    output zzz

Options:
    -h/--help
    -v/--version

Outputs:
    output zzz
'''

# # Parse options
def parse_options(argv):

    opts, _args = getopt.getopt(argv[1:], "hvVi:c:o:n:",
                                    ["help",
                                     "version",
                                     "input",
                                     "config",
                                     "output",
                                     "negative"])

    # Default working dir and config file
    input_file = ""
    output_folder = ""
    config_file = ""
    negative_file = None

    # Basic options
    for option, value in opts:
        if option in ("-h", "--help"):
            print(help_message)
            sys.exit(0)
        if option in ("-v", "-V", "--version"):
            print("xxx.py V%s" % (get_version()))
            sys.exit(0)
        if option in ("-i", "--input"):
            input_file = value
        if option in ("-o", "--output"):
            output_folder = value
        if option in ("-c", "--config"):
            config_file = value
        if option in ("-n", "--negative"):
            negative_file = value

    if input_file == "" or output_folder == "":
        print(help_message)
        sys.exit(0)

    output_folder = os.path.join(output_folder, '')

    return (input_file, config_file, output_folder, negative_file)

# # Decoy Reverse Forward protein
def protein_type(protein_sequence, lProtein=None):
    sProteins = protein_sequence.replace('{', '')
    sProteins = sProteins.replace('}', '')
    asProteins = sProteins.split(',')
    if lProtein != None:
        del lProtein[:]
        for sProtein in asProteins:
            sProtein = sProtein.strip()
            if sProtein not in lProtein:
                lProtein.extend(asProteins[:])
    for sProtein in asProteins:
        if not (sProtein.startswith(label_train_str) or sProtein.startswith(label_test_str)):
            return LabelFwd
    for sProtein in asProteins:
        if sProtein.startswith(label_test_str):
            return LabelTest
    if label_reserve_str != '':
        for sProtein in asProteins:
            if sProtein.startswith(label_reserve_str):
                return LabelRevReserve
    return LabelRevTrain

def protein_type_ecoli(protein_list):
    for one_protein in protein_list:
        if one_protein.startswith(label_ecoli_str):
            return LabelEcoli
    
    return LabelNotEcoli

# # split based on score agreement
def categorize_score_agreement(val):
    if val == 3:
        return '3'
    elif val == 2:
        return '2'
    else:
        return '1'

one_top_list = [0, 1, 2, 4]
# # split based on score agreement
def categorize_score_agreement_more_info(val):
    if val in one_top_list:
        return '1'
    elif val == 7:
        return 'A'
    elif val == 3:
        return 'XW'
    elif val == 5:
        return 'MW'
    elif val == 6:
        return 'MX'
    else:
        die("error")

two_top_list = [3, 5, 6]
# # split based on score agreement
def categorize_score_agreement_one_score_agreed(val):
    if val in two_top_list:
        return '2'
    elif val == 7:
        return '3'
    elif val == 0:
        return 'N'
    elif val == 1:
        return 'W'
    elif val == 2:
        return 'X'
    elif val == 4:
        return 'M'
    else:
        die("error")

# # split based on charge
def categorize_parent_charge(val):
    if val == 1:
        return '1'
    elif val == 2:
        return '2'
    elif val >= 3:
        return '3'

# # get the categories
def categorize(input_file, negative_file):

    psm_dict = {}
    # psm_set = set()
    
    # read line with csv
    f = open(input_file, 'r')
    
    # skip header
    _sHeader = f.readline()
    # get data
    for sLine in f:
        if sLine.startswith('#'):
            continue
        sLine = sLine.strip().split('\t')
        
        if len(sLine) == 32:
            PsmFields_obj = PsmFields2._make(sLine)
        elif len(sLine) == 35:
            PsmFields_obj = PsmFields4._make(sLine)
        else:
            PsmFields_obj = PsmFields3._make(sLine)
        sParentCharge = categorize_parent_charge(int(PsmFields_obj.ParentCharge))
        sScoreAgreement = categorize_score_agreement(int(PsmFields_obj.ScoreAgreement))
        # sScoreAgreement = categorize_score_agreement_more_info(int(PsmFields_obj.ScoreAgreement))
        # sScoreAgreement = categorize_score_agreement_one_score_agreed(int(PsmFields_obj.ScoreAgreement))

        sKey = sParentCharge + '_' + sScoreAgreement
        psm_obj = PSM(PsmFields_obj)
        '''
        if abs(psm_obj.fMassDiff) > mass_tolerance:
            continue
        '''
        '''
        if psm_obj.ScoreAgreement == 1 and psm_obj.iLocalRank > 0:
            continue
        '''
        '''
        if 'I' in PsmFields_obj.OriginalPeptide:
            psm_set.add(PsmFields_obj.OriginalPeptide.replace('I', 'L'))
        else:
            psm_set.add(PsmFields_obj.OriginalPeptide)
        '''
        if sKey not in psm_dict:
            psm_list = []
            psm_list.append(psm_obj)
            psm_dict[sKey] = psm_list
        else:
            psm_list = psm_dict[sKey]
            psm_list.append(psm_obj)
            
    f.close()
    
    psm_neg_list = []


    return (psm_dict, psm_neg_list)

# # Division error handling
divide = sipros_post_module.divide
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


# using the global rank
def get_cutoff_rank_product2(FDR_threshold, lPsm):
    F_list = []
    T_list = []
    TT_list = []
    S_list = []
    for oPsm in lPsm:
        if oPsm.RealLabel != LabelRevTrain:
            T_list.append(oPsm.fRankProduct)
            if oPsm.RealLabel == LabelTest:
                S_list.append(oPsm.fRankProduct)
            else:
                TT_list.append(oPsm.fRankProduct)
        else:
            F_list.append(oPsm.fRankProduct)
    F_list = np.array(F_list)
    T_list = np.array(T_list)
    TT_list = np.array(TT_list)
    S_list = np.array(S_list)
    prev_TF_num = 0
    final_cutoff_score = np.amin(T_list, axis=0)
    _final_accept = False
    T_num_best = 0
    TT_num_best = 0
    F_num_best = 0
    S_num_best = 0
    TF_num_best = 0
    for cutoff_score in T_list:
        F_num = (F_list <= cutoff_score).sum()
        T_num = (T_list <= cutoff_score).sum()
        TT_num = (TT_list <= cutoff_score).sum()
        S_num = (S_list <= cutoff_score).sum()
        # (F_num, T_num) = get_target_decoy_hit(lPsm, cutoff_score, iScoreIndex)
        TF_num = F_num + T_num
        TF_num = S_num + TT_num
        (FDR_accept, FDR_value) = FDR_calculator(S_num, TT_num)
        # update final values if conditions satisfies
        # 1) FDR_accept is True
        # 2) FDR_value should be less than or equal to FDR_threshold_1
        # 3) TF_num is greater than to previous TF_num
        if (FDR_accept is True) and (FDR_value <= FDR_threshold) and (TF_num > prev_TF_num) :
            final_cutoff_score = cutoff_score
            _final_accept = FDR_accept
            # previous TF_num
            prev_TF_num = TF_num
            T_num_best = T_num
            F_num_best = F_num
            S_num_best = S_num
            TF_num_best = TF_num

    print("Q-value_cutoff:\t%f\t%d\t%d\t%d" % (final_cutoff_score, T_num_best - S_num_best, F_num_best, S_num_best))
    return final_cutoff_score


# using the q value rank
def get_cutoff_q_rank_product(FDR_threshold, lPsm):
    iNumScores = lPsm[0].iNumScores

    T_num = 0
    F_num = 0
    S_num = 0
    for oPsm in lPsm:
        if oPsm.RealLabel != LabelRevTrain:
            T_num += 1
            if oPsm.RealLabel == LabelTest:
                S_num += 1
        else:
            F_num += 1
    print("Before Filtering:\t\t%d\t%d\t%d" % (T_num - S_num, F_num, S_num))

    for i in range(iNumScores):
        newlist = sorted(lPsm, key=lambda x: x.lfScores[i], reverse=True)
        T_num = 0
        F_num = 0
        for j in range(len(newlist)):
            if newlist[j].RealLabel != LabelRevTrain:
                T_num += 1
            else:
                F_num += 1
            (_FDR_accept, FDR_value) = FDR_calculator(F_num, T_num)
            newlist[j].lRanks.append(FDR_value)
        fSmallestQ = 1
        for j in range(len(newlist) - 1, -1, -1):
            if fSmallestQ > newlist[j].lRanks[i]:
                fSmallestQ = newlist[j].lRanks[i]
            if newlist[j].lRanks[i] > fSmallestQ:
                newlist[j].lRanks[i] = fSmallestQ

    for oPsm in lPsm:
        # fTemp = oPsm.lRanks[0] * oPsm.lRanks[1] * oPsm.lRanks[2]
        # oPsm.fRankProduct = np.power(float(fTemp), 1.0 / 3.0)
        fTemp = (1 - ((1 - oPsm.lRanks[0]) * (1 - oPsm.lRanks[1]) * (1 - oPsm.lRanks[2])))
        oPsm.fRankProduct = fTemp

    final_cutoff_score = get_cutoff_rank_product2(FDR_threshold, lPsm)
    return final_cutoff_score

def mark_training_label(lPsm, final_cutoff_score, count_list=None, sKey=None):
    num_positive = 0
    for oPsm in lPsm:
        if oPsm.RealLabel == LabelRevTrain:
            oPsm.TrainingLabel = LabelNegative
            '''
            if random.random() < 0.5:
                oPsm.TrainingLabel = LabelPositive
            '''
        else:
            '''
            if oPsm.ScoreAgreement < 2:
                if oPsm.iLocalRank == 0:
                    oPsm.TrainingLabel = LabelPositive
                else:
                    oPsm.TrainingLabel = LabelUnknown
                continue
            '''
            if oPsm.fRankProduct <= final_cutoff_score:
                oPsm.TrainingLabel = LabelPositive
                num_positive += 1
            else:
                oPsm.TrainingLabel = LabelUnknown
                oPsm.TrainingLabel = LabelPositive
                if sKey != None:
                    if sKey == '3_2' or sKey == '2_3' or sKey == '3_3':
                        oPsm.TrainingLabel = LabelUnknown
    '''            
    if num_positive >= 10000:
        list_sorted = sorted(lPsm, key=lambda x: x.fRankProduct)
        num_positive *= 0.95
        for oPsm in list_sorted:
            if num_positive > 0:
                if oPsm.TrainingLabel == LabelPositive:
                    num_positive -= 1
            else:
                if oPsm.TrainingLabel == LabelPositive:
                    oPsm.TrainingLabel = LabelUnknown
    '''
    if count_list != None:
        for oPsm in lPsm:
            if oPsm.fRankProduct <= final_cutoff_score:
                if oPsm.RealLabel == LabelFwd:
                    count_list[0] += 1
                elif oPsm.RealLabel == LabelRevTrain:
                    count_list[1] += 1
                else:
                    count_list[2] += 1

def show_TP_TN_FP_FN(label_np, predict_np):
    true_np = (label_np == predict_np)
    TP = label_np[true_np].sum()
    TN = (label_np[true_np] == 0).sum()
    false_np = (label_np != predict_np)
    FP = (label_np[false_np] == 0).sum()
    FN = label_np[false_np].sum()
    print("TP\tTN\tFP\tFN")
    print("%d\t%d\t%d\t%d" % (TP, TN, FP, FN))
    
def show_Fdr_category(psm_dict):
    Best_last_list = [0, 0, 0]
    for _key, lPsm in psm_dict.items():
        #print _key
        list_sorted = sorted(lPsm, key=lambda x: (x.fPredictProbability, 1 - x.fRankProduct) , reverse=True)
        T_num = 0
        F_num = 0
        Fwd_num = 0
        Rev_num = 0
        Shu_num = 0
        Best_list = [0, 0, 0]
        for oPsm in list_sorted:
            if oPsm.RealLabel == LabelFwd:
                T_num += 1
                Fwd_num += 1
            elif oPsm.RealLabel == LabelRevTrain:
                F_num += 1
                Rev_num += 1
            else:
                T_num += 1
                Shu_num += 1
            (FDR_accept, FDR_value) = FDR_calculator(Shu_num, Fwd_num)
            if (FDR_accept is True) and (FDR_value <= 0.01): #(Shu_num <= shuffle_num_control[sKey]):  # and (Shu_num <= 618) and (FDR_value <= 0.01)
                if Best_list[0] < Fwd_num:
                    Best_list = [Fwd_num, Rev_num, Shu_num]
        #sys.stdout.write('%d\t%d\t%d\n' % (Best_list[0], Best_list[1], Best_list[2]))
        Best_last_list = [Best_last_list[0] + Best_list[0], Best_last_list[1] + Best_list[1], Best_last_list[2] + Best_list[2]]
    sys.stdout.write('%d\t%d\t%d\t' % (Best_last_list[0], Best_last_list[1], Best_last_list[2]))
    pass

def filter_Fdr(psm_list, fdr):
    list_sorted = sorted(psm_list, key=lambda x: (x.fPredictProbability, 1 - x.fRankProduct) , reverse=True)
    T_num = 0
    F_num = 0
    Fwd_num = 0
    Rev_num = 0
    Shu_num = 0
    Best_list = [0, 0, 0]
    psm_filtered_list = []
    psm_left_list = []
    cutoff_probability = 0.0
    for oPsm in list_sorted:
        if oPsm.TrainingLabel == LabelFiltered:
            continue
        if oPsm.RealLabel == LabelFwd:
            T_num += 1
            Fwd_num += 1
        elif oPsm.RealLabel == LabelRevTrain:
            F_num += 1
            Rev_num += 1
        else:
            T_num += 1
            Shu_num += 1
        (FDR_accept, FDR_value) = FDR_calculator(Shu_num, Fwd_num)
        if (FDR_accept is True) and (FDR_value <= fdr): #(Shu_num <= shuffle_num_control[sKey]):  # and (Shu_num <= 618) and (FDR_value <= 0.01)
            if (Best_list[0] + Best_list[2]) < (Fwd_num + Shu_num):
                Best_list = [Fwd_num, Rev_num, Shu_num]
                cutoff_probability = oPsm.fPredictProbability
    
    sys.stdout.write('\n%d\t%d\t%d\t\n' % (Best_list[0], Best_list[1], Best_list[2]))
    
    for oPsm in list_sorted:
        if oPsm.TrainingLabel == LabelFiltered:
            continue
        if oPsm.fPredictProbability >= cutoff_probability:
            psm_filtered_list.append(oPsm)
            oPsm.TrainingLabel = LabelFiltered
        else:
            psm_left_list.append(oPsm)
    
    return (psm_filtered_list, psm_left_list)

def filter_Fdr2(psm_list, fdr):
    list_sorted = sorted(psm_list, key=lambda x: (x.fPredictProbability, 1 - x.fRankProduct) , reverse=True)
    T_num = 0
    F_num = 0
    Fwd_num = 0
    Rev_num = 0
    Shu_num = 0
    Best_list = [0, 0, 0]
    psm_filtered_list = []
    psm_left_list = []
    cutoff_probability = 0.0
    for oPsm in list_sorted:
        if oPsm.TrainingLabel == LabelFiltered:
            continue
        if oPsm.fRtPvalue < 0.0001:
            continue
        if oPsm.RealLabel == LabelFwd:
            T_num += 1
            Fwd_num += 1
        elif oPsm.RealLabel == LabelRevTrain:
            F_num += 1
            Rev_num += 1
        else:
            T_num += 1
            Shu_num += 1
        (FDR_accept, FDR_value) = FDR_calculator(Shu_num, Fwd_num)
        if (FDR_accept is True) and (FDR_value <= fdr): #(Shu_num <= shuffle_num_control[sKey]):  # and (Shu_num <= 618) and (FDR_value <= 0.01)
            if (Best_list[0] + Best_list[2]) < (Fwd_num + Shu_num):
                Best_list = [Fwd_num, Rev_num, Shu_num]
                cutoff_probability = oPsm.fPredictProbability
    
    sys.stdout.write('\n%d\t%d\t%d\t\n' % (Best_list[0], Best_list[1], Best_list[2]))
    
    for oPsm in list_sorted:
        if oPsm.TrainingLabel == LabelFiltered:
            continue
        if oPsm.fPredictProbability >= cutoff_probability:
            psm_filtered_list.append(oPsm)
            oPsm.TrainingLabel = LabelFiltered
        else:
            psm_left_list.append(oPsm)
    
    return (psm_filtered_list, psm_left_list)

def show_Fdr_charge(psm_list):
    psm_new_list_1 = []
    psm_new_list_2 = []
    for oPsm in psm_list:
        if oPsm.ParentCharge <= 2:
            psm_new_list_1.append(oPsm)
        else:
            psm_new_list_2.append(oPsm)
    show_Fdr(psm_new_list_1, None, None)
    show_Fdr(psm_new_list_2, None, None)

psm_num_list = []
true_psm_num_list = []
pep_num_list = []
true_pep_num_list = []

def show_Fdr_varied(psm_list, fdr):
    
    # list_sorted = sorted(psm_list, key=lambda x: (x.fPredictProbability, 1 - x.fRankProduct) , reverse=True)
    list_sorted = sorted(psm_list, key=lambda x: (x.fPredictProbability) , reverse=True)
    T_num = 0
    F_num = 0
    Fwd_num = 0
    train_num = 0
    rev_num = 0
    ecoli_num = 0
    best_ecoli_num = 0
    Best_list = [0, 0, 0]

    psm_filtered_list = []
    
    cutoff_probability = 0.0
    
    peptide_set = set()
    num_fwr_pep = 0
    num_rev_pep = 0
    num_ecoli_pep = 0
    best_ecoli_pep = 0
    best_fwr_pep = 0
    best_rev_pep = 0
    
    # without considering training label
    for oPsm in list_sorted:
        if oPsm.TrainingLabel == LabelFiltered:
            continue
        
        pep_str = oPsm.IdentifiedPeptide + '_' + str(oPsm.ParentCharge)
        if pep_str not in peptide_set:
            if oPsm.RealLabel == LabelFwd:
                num_fwr_pep += 1
                peptide_set.add(pep_str)
                if oPsm.ecoli_label == LabelEcoli:
                    num_ecoli_pep += 1
            elif oPsm.RealLabel == LabelTest:
                num_rev_pep += 1
                peptide_set.add(pep_str)
        
        if oPsm.RealLabel == LabelFwd:
            T_num += 1
            Fwd_num += 1
            if oPsm.ecoli_label == LabelEcoli:
                ecoli_num += 1
        elif oPsm.RealLabel == LabelRevTrain:
            F_num += 1
            train_num += 1
        elif oPsm.RealLabel == LabelTest:
            T_num += 1
            rev_num += 1
        (FDR_accept, FDR_value) = FDR_calculator(rev_num, Fwd_num)
        if (FDR_accept is True) and (FDR_value <= fdr): #(rev_num <= shuffle_num_control[sKey]):  # and (rev_num <= 618) and (FDR_value <= 0.01)
            if (Best_list[0] + Best_list[2]) < (Fwd_num + rev_num):
                Best_list = [Fwd_num, train_num, rev_num]
                best_ecoli_num = ecoli_num
                cutoff_probability = oPsm.fPredictProbability
        (FDR_accept, FDR_value) = FDR_calculator(num_rev_pep, num_fwr_pep)
        if (FDR_accept is True) and (FDR_value <= fdr): #(rev_num <= shuffle_num_control[sKey]):  # and (rev_num <= 618) and (FDR_value <= 0.01)
            if (best_fwr_pep + best_rev_pep) < (num_fwr_pep + num_rev_pep):        
                best_fwr_pep = num_fwr_pep
                best_rev_pep = num_rev_pep
                best_ecoli_pep = num_ecoli_pep
            
    for oPsm in list_sorted:
        if oPsm.fPredictProbability >= cutoff_probability:
            psm_filtered_list.append(oPsm)
    
    # sys.stdout.write('\n'+str(len(psm_filtered_list))+'\t')
    
    print("{:,d} ({:,d}) ({:.2f}%)\t{:,d} ({:,d}) ({:.2f}%)".format(Best_list[0], best_ecoli_num, (100* float(Best_list[2])/float(Best_list[0])), best_fwr_pep, best_ecoli_pep, (100.0*float(best_rev_pep)/float(best_fwr_pep))))
    psm_num_list.append(Best_list[0])
    true_psm_num_list.append(best_ecoli_num)
    pep_num_list.append(best_fwr_pep)
    true_pep_num_list.append(best_ecoli_pep)
    return psm_filtered_list

def show_Fdr(psm_list, sKey, fdr=None):
    
    fdr_float = 0.01
    if fdr !=None:
        fdr_float = fdr
    
    # list_sorted = sorted(psm_list, key=lambda x: (x.fPredictProbability, 1 - x.fRankProduct) , reverse=True)
    list_sorted = sorted(psm_list, key=lambda x: (x.fPredictProbability) , reverse=True)
    Fwd_num = 0
    Rev_num = 0
    Shu_num = 0
    Best_list = [0, 0, 0]
    '''
    fwr = open("fwr.txt", "a")
    rev = open("rev.txt", "a")
    '''
    
    peptide_set = set()
    num_fwr_pep = 0
    num_shu_pep = 0
    best_fwr_pep = 0
    best_shu_pep = 0
    
    psm_filtered_list = []
    cutoff_probability = 0.0
    # without considering training label
    for oPsm in list_sorted:
        if oPsm.TrainingLabel == LabelFiltered:
            continue
        
        pep_str = oPsm.IdentifiedPeptide + '_' + str(oPsm.ParentCharge)
        if pep_str not in peptide_set:
            if oPsm.RealLabel == LabelFwd:
                num_fwr_pep += 1
                peptide_set.add(pep_str)
            elif oPsm.RealLabel == LabelTest:
                num_shu_pep += 1
                peptide_set.add(pep_str)
        
        if oPsm.RealLabel == LabelFwd:
            Fwd_num += 1
        elif oPsm.RealLabel == LabelRevTrain:
            Rev_num += 1
        elif oPsm.RealLabel == LabelTest:
            Shu_num += 1
        else:
            print('error 768.')
        (FDR_accept, FDR_value) = FDR_calculator(Shu_num, Fwd_num)
        if (FDR_accept is True) and (FDR_value <= fdr_float): #(Shu_num <= shuffle_num_control[sKey]):  # and (Shu_num <= 618) and (FDR_value <= 0.01)
            if (Best_list[0] + Best_list[2]) < (Fwd_num + Shu_num):
                Best_list = [Fwd_num, Rev_num, Shu_num]
                cutoff_probability = oPsm.fPredictProbability
                best_fwr_pep = num_fwr_pep
                best_shu_pep = num_shu_pep
            

    for oPsm in list_sorted:
        if oPsm.fPredictProbability >= cutoff_probability:
            psm_filtered_list.append(oPsm)
     
    sys.stdout.write('\t'+str(len(psm_filtered_list))+'\n')
    sys.stdout.write("{:,d}\n{:.2f}%\n{:.2f}%\n{:,d}\n{:.2f}%\n".format(Best_list[0], 100.0*float(Best_list[1])/float(Best_list[0]), 100.0*float(Best_list[2])/float(Best_list[0]), best_fwr_pep, 100.0*float(best_shu_pep)/float(best_fwr_pep)))
    # sys.stdout.write("%d\t[%.2f%%]\t[%.2f%%]\t%d\t[%.2f%%]\n" % (Best_list[0], 100.0*float(Best_list[2])/float(Best_list[0]), 100.0*float(Best_list[1])/float(Best_list[0]), best_fwr_pep, 100.0*float(best_shu_pep)/float(best_fwr_pep)))
    
    return psm_filtered_list

def show_Fdr_Pep(psm_list, sKey, fdr=None):
    
    fdr_float = 0.01
    if fdr !=None:
        fdr_float = fdr
    
    # list_sorted = sorted(psm_list, key=lambda x: (x.fPredictProbability, 1 - x.fRankProduct) , reverse=True)
    list_sorted = sorted(psm_list, key=lambda x: (x.fPredictProbability) , reverse=True)
    T_num = 0
    F_num = 0
    Fwd_num = 0
    Rev_num = 0
    Shu_num = 0
    Best_list = [0, 0, 0]
    '''
    fwr = open("fwr.txt", "a")
    rev = open("rev.txt", "a")
    '''
    
    peptide_set = set()
    num_fwr_pep = 0
    num_shu_pep = 0
    best_fwr_pep = 0
    best_shu_pep = 0
    
    psm_filtered_list = []
    cutoff_probability = 0.0
    # without considering training label
    for oPsm in list_sorted:
        if oPsm.TrainingLabel == LabelFiltered:
            continue
        
        pep_str = oPsm.IdentifiedPeptide + '_' + str(oPsm.ParentCharge)
        if pep_str not in peptide_set:
            if oPsm.RealLabel == LabelFwd:
                num_fwr_pep += 1
                peptide_set.add(pep_str)
            elif oPsm.RealLabel == LabelTest:
                num_shu_pep += 1
                peptide_set.add(pep_str)
        
        if oPsm.RealLabel == LabelFwd:
            T_num += 1
            Fwd_num += 1
        elif oPsm.RealLabel == LabelRevTrain:
            F_num += 1
            Rev_num += 1
        elif oPsm.RealLabel == LabelTest:
            T_num += 1
            Shu_num += 1
        (FDR_accept, FDR_value) = FDR_calculator(num_shu_pep, num_fwr_pep)
        if (FDR_accept is True) and (FDR_value <= fdr_float): #(Shu_num <= shuffle_num_control[sKey]):  # and (Shu_num <= 618) and (FDR_value <= 0.01)
            if (Best_list[0] + Best_list[2]) < (Fwd_num + Shu_num):
                Best_list = [Fwd_num, Rev_num, Shu_num]
                cutoff_probability = oPsm.fPredictProbability
                best_fwr_pep = num_fwr_pep
                best_shu_pep = num_shu_pep
            

    for oPsm in list_sorted:
        if oPsm.fPredictProbability >= cutoff_probability:
            psm_filtered_list.append(oPsm)
     
    sys.stdout.write('\t'+str(len(psm_filtered_list))+'\n')
    
    sys.stdout.write("%d\t[%.2f%%]\t[%.2f%%]\t%d\t[%.2f%%]\n" % (Best_list[0], 100.0*float(Best_list[2])/float(Best_list[0]), 100.0*float(Best_list[1])/float(Best_list[0]), best_fwr_pep, 100.0*float(best_shu_pep)/float(best_fwr_pep)))
    
    return psm_filtered_list

def show_Fdr_peptide(psm_list, sKey, fdr=None):
    
    fdr_float = 0.01
    if fdr !=None:
        fdr_float = fdr
    
    # list_sorted = sorted(psm_list, key=lambda x: (x.fPredictProbability, 1 - x.fRankProduct) , reverse=True)
    list_sorted = sorted(psm_list, key=lambda x: (x.fPredictProbability) , reverse=True)
    T_num = 0
    F_num = 0
    Fwd_num = 0
    Rev_num = 0
    Shu_num = 0
    Best_list = [0, 0, 0]

    psm_filtered_list = []
    cutoff_probability = 0.0
    peptide_set = set()
    # without considering training label
    for oPsm in list_sorted:
        pep_str = oPsm.IdentifiedPeptide + '_' + str(oPsm.ParentCharge)
        if pep_str in peptide_set:
            continue
        else:
            peptide_set.add(pep_str)
        if oPsm.TrainingLabel == LabelFiltered:
            continue
        if oPsm.RealLabel == LabelFwd:
            T_num += 1
            Fwd_num += 1
        elif oPsm.RealLabel == LabelRevTrain:
            F_num += 1
            Rev_num += 1
        else:
            T_num += 1
            Shu_num += 1
        (FDR_accept, FDR_value) = FDR_calculator(Shu_num, Fwd_num)
        if (FDR_accept is True) and (FDR_value <= fdr_float): #(Shu_num <= shuffle_num_control[sKey]):  # and (Shu_num <= 618) and (FDR_value <= 0.01)
            if (Best_list[0] + Best_list[2]) < (Fwd_num + Shu_num):
                Best_list = [Fwd_num, Rev_num, Shu_num]
                cutoff_probability = oPsm.fPredictProbability
            

    for oPsm in list_sorted:
        if oPsm.fPredictProbability >= cutoff_probability:
            psm_filtered_list.append(oPsm)
   
    sys.stdout.write('\t'+str(len(psm_filtered_list))+'\n')
    
    sys.stdout.write('%d\t%d\t%d\t' % (Best_list[0], Best_list[1], Best_list[2]))
    
    return psm_filtered_list

def get_RT_data(psm_filtered_list, output_folder, psm_list):
    
    # create a temporary directory
    if not os.path.exists(output_folder + 'temp'):
        os.makedirs(output_folder + 'temp')
    else:
        for the_file in os.listdir(output_folder + 'temp'):
            file_path = os.path.join(output_folder + 'temp', the_file)
            try:
                if os.path.isfile(file_path):
                    os.unlink(file_path)
                # elif os.path.isdir(file_path): shutil.rmtree(file_path)
            except Exception as e:
                print(e)
    output_folder += 'temp/'
    
    shift_threshold_float = 7.5
    
    # shift the RT
    ft2_filename_list = []
    ft2_filename_max_rt_dict = {}
    for oPsm in psm_list:
        if oPsm.FileName not in ft2_filename_max_rt_dict:
            ft2_filename_list.append(oPsm.FileName)
            ft2_filename_max_rt_dict[oPsm.FileName] = oPsm.fRtMeasured
        else:
            if oPsm.fRtMeasured > ft2_filename_max_rt_dict[oPsm.FileName]:
                ft2_filename_max_rt_dict[oPsm.FileName] = oPsm.fRtMeasured
    
    ft2_filename_list = sorted(ft2_filename_list)
    
    for oPsm in psm_list:
        if oPsm.fRtMeasured <= shift_threshold_float and oPsm.FileName != ft2_filename_list[0]:
            oPsm.bFileNameChanged = True
            ind_filename_int = ft2_filename_list.index(oPsm.FileName)
            oPsm.FileName = ft2_filename_list[ind_filename_int - 1]
            oPsm.fRtMeasured += ft2_filename_max_rt_dict[oPsm.FileName]
            oPsm.sRTime = "%.5f" % oPsm.fRtMeasured
    
    # Prepare the training data
    # # write into hard disk
    psm_filtered_list = sorted(psm_filtered_list, key=lambda x: (x.fPredictProbability) , reverse=True)
    ft2_train_test_model_result_list = []
    train_file_dict = {}
    test_file_dict = {}
    train_num_dict = {}
    for oPsm in psm_filtered_list:
        if oPsm.RealLabel == LabelRevTrain or oPsm.RealLabel == LabelTest:
            continue
        if oPsm.FileName in train_file_dict:
            if train_num_dict[oPsm.FileName] < 1000:
                fw = train_file_dict[oPsm.FileName]
                fw.write((oPsm.IdentifiedPeptide[1:-1]).replace('~', '(Oxidation)'))
                fw.write('\t')
                fw.write(oPsm.sRTime)
                fw.write('\n')
                train_num_dict[oPsm.FileName] += 1
        else:
            train_file_str = output_folder + oPsm.FileName.split('.')[0] + '_train.txt'
            test_file_str = output_folder + oPsm.FileName.split('.')[0] + '_test.txt'
            model_file_str = output_folder + oPsm.FileName.split('.')[0] + '_model.txt'
            result_file_str = output_folder + oPsm.FileName.split('.')[0] + '_result.csv'
            ft2_train_test_model_result_list.append((oPsm.FileName, train_file_str, test_file_str, model_file_str, result_file_str))
            fw_test = open(test_file_str, 'w')
            test_file_dict[oPsm.FileName] = fw_test
            fw = open(train_file_str, 'w')
            train_file_dict[oPsm.FileName] = (fw)
            fw.write((oPsm.IdentifiedPeptide[1:-1]).replace('~', '(Oxidation)'))
            fw.write('\t')
            fw.write(oPsm.sRTime)
            fw.write('\n')
            train_num_dict[oPsm.FileName] = 1
    
    for _key, value in train_file_dict.items():
        value.close()
    
    
    # for oPsm in psm_list:
    for oPsm in psm_filtered_list:
        if oPsm.FileName in test_file_dict:
            fw_test = test_file_dict[oPsm.FileName]
            fw_test.write((oPsm.IdentifiedPeptide[1:-1]).replace('~', '(Oxidation)'))
            fw_test.write('\n')
    
               
    for _key, value in test_file_dict.items():
        value.close()
        
    PGM = "/home/xgo/ORNL/Sipros/MachineLearning/RetentionTime/OpenMS/OpenMS-build/bin/"
    iNumThreads = cpu_count() - 1
    qUnprocessed = Queue(len(ft2_train_test_model_result_list)*2)
    qProcessed = Queue(len(ft2_train_test_model_result_list)*2)
    
    # run RTModel
    for (_ft_file_str, train_file_str, test_file_str, model_file_str, result_file_str) in ft2_train_test_model_result_list:
        call_list = [PGM+"RTModel", '-in', train_file_str, '-out', model_file_str, '-cv:skip_cv']
        qUnprocessed.put(call_list)
    for _i in range(iNumThreads):
        PsmProcessor = RTModel(qUnprocessed, qProcessed)
        PsmProcessor.daemon = True
        PsmProcessor.start()
        qUnprocessed.put(None)
    
    iNumRankers = iNumThreads
    while True:
        call_list = qProcessed.get(True)
        if call_list is None:
            iNumRankers -= 1
            if iNumRankers == 0:
                break
            else:
                continue
    '''
    '''
    # run RTPredict
    for (_ft_file_str, train_file_str, test_file_str, model_file_str, result_file_str) in ft2_train_test_model_result_list:
        call_list = [PGM+"RTPredict", '-in_text', test_file_str, '-svm_model', model_file_str, '-out_text:file', result_file_str]
        qUnprocessed.put(call_list)
    for _i in range(iNumThreads):
        PsmProcessor = RTPredict(qUnprocessed, qProcessed)
        PsmProcessor.daemon = True
        PsmProcessor.start()
        qUnprocessed.put(None)
    
    iNumRankers = iNumThreads
    while True:
        call_list = qProcessed.get(True)
        if call_list is None:
            iNumRankers -= 1
            if iNumRankers == 0:
                break
            else:
                continue
    
    # collect RTtime
    ft_pep_rtime_dict = {}
    for (ft_file_str, train_file_str, test_file_str, model_file_str, result_file_str) in ft2_train_test_model_result_list:
        pep_rtime_dict = {}
        ft_pep_rtime_dict[ft_file_str] = pep_rtime_dict
        with open(result_file_str, 'r') as fr:
            for line_str in fr:
                words = line_str.strip().split()
                pep_str = '[' + words[0].replace('(Oxidation)', '~') + ']'
                if pep_str in pep_rtime_dict:
                    print('error in pep_str in pep_rtime_dict')
                else:
                    rt_float = float(words[1])
                    pep_rtime_dict[pep_str] = rt_float
    
    # assign rtime to psm, and roll back the file name
    psm_na_rt_list = []
    # for oPsm in psm_list:
    for oPsm in psm_filtered_list:
        if oPsm.FileName not in ft_pep_rtime_dict:
            psm_na_rt_list.append(oPsm)
        else:
            pep_rtime_dict = ft_pep_rtime_dict[oPsm.FileName]
            oPsm.fRtPredict = float(pep_rtime_dict[oPsm.IdentifiedPeptide])
        
        if oPsm.bFileNameChanged:
            idx_filename_int = ft2_filename_list.index(oPsm.FileName)
            oPsm.FileName = ft2_filename_list[idx_filename_int + 1]
            oPsm.bFileNameChanged = False
            
    # save the measure time and predict time and label
    out_folder_str = "/media/xgo/Seagate/Proteomics/Experiments/BenchmarkRev/Sipros10/0820/fig3Sumplement/"
    out_fw_dict = {}
    # for oPsm in psm_list:
    for oPsm in psm_filtered_list:
        if oPsm.FileName in out_fw_dict:
            fw = out_fw_dict[oPsm.FileName]
            fw.write(str(oPsm.RealLabel))
            fw.write("\t")
            fw.write(str(oPsm.sRTime))
            fw.write("\t")
            fw.write(str(oPsm.fRtPredict))
            fw.write('\n')
        else:
            fw = open(out_folder_str + oPsm.FileName.split('.')[0] + '.txt', 'w')
            fw.write("RealLabel\tMeasuredTime\tPredictedTime\n")
            fw.write(str(oPsm.RealLabel))
            fw.write("\t")
            fw.write(str(oPsm.sRTime))
            fw.write("\t")
            fw.write(str(oPsm.fRtPredict))
            fw.write('\n')
            out_fw_dict[oPsm.FileName] = fw
    for _filename, fw in out_fw_dict.items():
        fw.close()
    return
        
    # filtering based on RT time
    range_list = []
    threshold_list = []
    distribuation_list = []
    num_range = 6
    scale_range = 20
    for i in range(1, num_range+1):
        range_list.append([])
        threshold_list.append(i* scale_range)
    threshold_list[num_range - 1] = 125
    for oPsm in psm_filtered_list:
        idx = 0
        for i in range(num_range):
            if oPsm.fRtMeasured <= threshold_list[i]:
                idx = i
                break
        range_list[idx].append(abs(oPsm.fRtMeasured - oPsm.fRtPredict))
        
    for value_list in range_list:
        value_list_np = np.array(value_list)
        std_float = np.std(value_list_np)
        mean_float = np.mean(value_list_np)
        distribuation_list.append((mean_float, std_float))
    
    for oPsm in psm_list:
        idx = 0
        for i in range(num_range):
            if oPsm.fRtMeasured <= threshold_list[i]:
                idx = i
                break
        oPsm.fRtPvalue = spst.norm(distribuation_list[idx][0], distribuation_list[idx][1]).cdf(abs(oPsm.fRtMeasured - oPsm.fRtPredict))
        oPsm.fRtPvalue = 1 - oPsm.fRtPvalue
    

def get_RT_elude_data(psm_filtered_list, output_folder, psm_list):

    # create a temporary directory
    if not os.path.exists(output_folder + 'temp'):
        os.makedirs(output_folder + 'temp')
    else:
        for the_file in os.listdir(output_folder + 'temp'):
            file_path = os.path.join(output_folder + 'temp', the_file)
            try:
                if os.path.isfile(file_path):
                    os.unlink(file_path)
                # elif os.path.isdir(file_path): shutil.rmtree(file_path)
            except Exception as e:
                print(e)
    output_folder += 'temp/'

    shift_threshold_float = 7.5
    
    # shift the RT
    ft2_filename_list = []
    ft2_filename_max_rt_dict = {}
    for oPsm in psm_list:
        if oPsm.FileName not in ft2_filename_max_rt_dict:
            ft2_filename_list.append(oPsm.FileName)
            ft2_filename_max_rt_dict[oPsm.FileName] = oPsm.fRtMeasured
        else:
            if oPsm.fRtMeasured > ft2_filename_max_rt_dict[oPsm.FileName]:
                ft2_filename_max_rt_dict[oPsm.FileName] = oPsm.fRtMeasured
    
    ft2_filename_list = sorted(ft2_filename_list)
    
    for oPsm in psm_list:
        if oPsm.fRtMeasured <= shift_threshold_float and oPsm.FileName != ft2_filename_list[0]:
            oPsm.bFileNameChanged = True
            ind_filename_int = ft2_filename_list.index(oPsm.FileName)
            oPsm.FileName = ft2_filename_list[ind_filename_int - 1]
            oPsm.fRtMeasured += ft2_filename_max_rt_dict[oPsm.FileName]
            oPsm.sRTime = "%.5f" % oPsm.fRtMeasured
    
    psm_filtered_list = sorted(psm_filtered_list, key=lambda x: (x.fPredictProbability) , reverse=True)
    ft2_train_test_model_result_list = []
    train_file_dict = {}
    test_file_dict = {}
    train_num_dict = {}
    for oPsm in psm_filtered_list:
        if oPsm.RealLabel == LabelRevTrain or oPsm.RealLabel == LabelTest:
            continue
        if oPsm.FileName in train_file_dict:
            if train_num_dict[oPsm.FileName] <= 1000:
                fw = train_file_dict[oPsm.FileName]
                fw.write((oPsm.IdentifiedPeptide[1:-1]).replace('~', '[Oxidation]'))
                fw.write('\t')
                fw.write(oPsm.sRTime)
                fw.write('\n')
                train_num_dict[oPsm.FileName] += 1
        else:
            train_file_str = output_folder + oPsm.FileName.split('.')[0] + '_train.txt'
            test_file_str = output_folder + oPsm.FileName.split('.')[0] + '_test.txt'
            model_file_str = output_folder + oPsm.FileName.split('.')[0] + '_model.txt'
            result_file_str = output_folder + oPsm.FileName.split('.')[0] + '_result.txt'
            ft2_train_test_model_result_list.append((oPsm.FileName, train_file_str, test_file_str, model_file_str, result_file_str))
            fw_test = open(test_file_str, 'w')
            test_file_dict[oPsm.FileName] = fw_test
            fw = open(train_file_str, 'w')
            train_file_dict[oPsm.FileName] = (fw)
            fw.write((oPsm.IdentifiedPeptide[1:-1]).replace('~', '[Oxidation]'))
            fw.write('\t')
            fw.write(oPsm.sRTime)
            fw.write('\n')
            train_num_dict[oPsm.FileName] = 1
    
    for _key, value in train_file_dict.items():
        value.close()
    
    # for oPsm in psm_list:
    for oPsm in psm_filtered_list:
        if oPsm.FileName in test_file_dict:
            fw_test = test_file_dict[oPsm.FileName]
            fw_test.write((oPsm.IdentifiedPeptide[1:-1]).replace('~', '[Oxidation]'))
            fw_test.write('\n')
            
    for _key, value in test_file_dict.items():
        value.close()
        
    PGM = "/usr/bin/"
    iNumThreads = cpu_count() - 1
    qUnprocessed = Queue(len(ft2_train_test_model_result_list)*2)
    qProcessed = Queue(len(ft2_train_test_model_result_list)*2)
    
    # run Elude
    for (_ft_file_str, train_file_str, test_file_str, model_file_str, result_file_str) in ft2_train_test_model_result_list:
        call_list = [PGM+"elude", '-t', train_file_str, '-s', model_file_str]
        qUnprocessed.put(call_list)
    for _i in range(iNumThreads):
        PsmProcessor = RTModel(qUnprocessed, qProcessed)
        PsmProcessor.daemon = True
        PsmProcessor.start()
        qUnprocessed.put(None)
    
    iNumRankers = iNumThreads
    while True:
        call_list = qProcessed.get(True)
        if call_list is None:
            iNumRankers -= 1
            if iNumRankers == 0:
                break
            else:
                continue
    '''
    '''
    # run RTPredict
    for (_ft_file_str, train_file_str, test_file_str, model_file_str, result_file_str) in ft2_train_test_model_result_list:
        call_list = [PGM+"elude", '-e', test_file_str, '-l', model_file_str, '-o', result_file_str]
        qUnprocessed.put(call_list)
    for _i in range(iNumThreads):
        PsmProcessor = RTPredict(qUnprocessed, qProcessed)
        PsmProcessor.daemon = True
        PsmProcessor.start()
        qUnprocessed.put(None)
    
    iNumRankers = iNumThreads
    while True:
        call_list = qProcessed.get(True)
        if call_list is None:
            iNumRankers -= 1
            if iNumRankers == 0:
                break
            else:
                continue
    
    # collect RTtime
    ft_pep_rtime_dict = {}
    for (ft_file_str, train_file_str, test_file_str, model_file_str, result_file_str) in ft2_train_test_model_result_list:
        pep_rtime_dict = {}
        ft_pep_rtime_dict[ft_file_str] = pep_rtime_dict
        with open(result_file_str, 'r') as fr:
            fr.next()
            fr.next()
            fr.next()
            for line_str in fr:
                words = line_str.strip().split()
                pep_str = '[' + words[0].replace('[Oxidation]', '~') + ']'
                if pep_str in pep_rtime_dict:
                    pass
                    # print 'error in pep_str in pep_rtime_dict'
                else:
                    rt_float = float(words[1])
                    pep_rtime_dict[pep_str] = rt_float
    
    # assign rtime to psm, and roll back the file name
    psm_na_rt_list = []
    for oPsm in psm_filtered_list:
        if oPsm.FileName not in ft_pep_rtime_dict:
            psm_na_rt_list.append(oPsm)
        else:
            pep_rtime_dict = ft_pep_rtime_dict[oPsm.FileName]
            if oPsm.IdentifiedPeptide not in pep_rtime_dict:
                print('error in pep_rtime_dict.')
            oPsm.fRtPredict = float(pep_rtime_dict[oPsm.IdentifiedPeptide])
        
        if oPsm.bFileNameChanged:
            idx_filename_int = ft2_filename_list.index(oPsm.FileName)
            oPsm.FileName = ft2_filename_list[idx_filename_int + 1]
            oPsm.bFileNameChanged = False
    
    if len(psm_na_rt_list) > 0:
        print('error in psm_na_rt_list.')
    
    # save the measure time and predict time and label
    out_folder_str = "/media/xgo/Seagate/Proteomics/Experiments/BenchmarkRev/Sipros10/0820/fig3Sumplement/"
    out_fw_dict = {}
    for oPsm in psm_filtered_list:
        if oPsm.FileName in out_fw_dict:
            fw = out_fw_dict[oPsm.FileName]
            fw.write(str(oPsm.RealLabel))
            fw.write("\t")
            fw.write(str(oPsm.sRTime))
            fw.write("\t")
            fw.write(str(oPsm.fRtPredict))
            fw.write('\n')
        else:
            fw = open(out_folder_str + oPsm.FileName.split('.')[0] + '.txt', 'w')
            fw.write("RealLabel\tMeasuredTime\tPredictedTime\n")
            fw.write(str(oPsm.RealLabel))
            fw.write("\t")
            fw.write(str(oPsm.sRTime))
            fw.write("\t")
            fw.write(str(oPsm.fRtPredict))
            fw.write('\n')
            out_fw_dict[oPsm.FileName] = fw
    for _filename, fw in out_fw_dict.items():
        fw.close()
    return    
    
    # filtering based on RT time
    range_list = []
    threshold_list = []
    distribuation_list = []
    num_range = 6
    scale_range = 20
    for i in range(1, num_range+1):
        range_list.append([])
        threshold_list.append(i* scale_range)
    threshold_list[num_range - 1] = 125
    for oPsm in psm_filtered_list:
        idx = 0
        for i in range(num_range):
            if oPsm.fRtMeasured <= threshold_list[i]:
                idx = i
                break
        range_list[idx].append(abs(oPsm.fRtMeasured - oPsm.fRtPredict))
        
    for value_list in range_list:
        value_list_np = np.array(value_list)
        std_float = np.std(value_list_np)
        mean_float = np.mean(value_list_np)
        distribuation_list.append((mean_float, std_float))
    
    for oPsm in psm_list:
        idx = 0
        for i in range(num_range):
            if oPsm.fRtMeasured <= threshold_list[i]:
                idx = i
                break
        oPsm.fRtPvalue = spst.norm(distribuation_list[idx][0], distribuation_list[idx][1]).cdf(abs(oPsm.fRtMeasured - oPsm.fRtPredict))
        oPsm.fRtPvalue = 1 - oPsm.fRtPvalue

# # thread class for ranking the PSM
class RTModel(Process):

    def __init__(self, qUnprocessed, qProcessed):
        super(RTModel, self).__init__()
        self.qUnprocessed = qUnprocessed
        self.qProcessed = qProcessed
        return

    def run(self):
        while True:
            call_list = self.qUnprocessed.get(True)
            if call_list is None:
                break
            call(call_list)
        self.qProcessed.put(None)
        return

# # thread class for ranking the PSM
class RTPredict(Process):

    def __init__(self, qUnprocessed, qProcessed):
        super(RTPredict, self).__init__()
        self.qUnprocessed = qUnprocessed
        self.qProcessed = qProcessed
        return

    def run(self):
        while True:
            call_list = self.qUnprocessed.get(True)
            if call_list is None:
                break
            call(call_list)
        self.qProcessed.put(None)
        return


# from keras.models import Sequential    
# from keras.layers import Dense, Dropout

def test_DeepLearning(psm_dict, psm_list, psm_neg_list, fdr_given=None, psm_left_list=None):
    # machine learning
    # # construct training data
    train_data_list = []
    train_label_list = []
    test_data_list = []
    positive_int = 1
    negative_int = 0
    psm_rank_list = []     
    for oPsm in psm_list:
        if len(oPsm.feature_list) == 0:
            print('check')
        if oPsm.TrainingLabel == LabelFiltered:
            continue
        if oPsm.RealLabel == LabelRevReserve:
            # train_data_list.append(oPsm.feature_list)
            # train_label_list.append(negative_int)
            continue
        test_data_list.append(oPsm.feature_list)
        psm_rank_list.append(oPsm)     
        if oPsm.RealLabel != LabelRevTrain and (oPsm.iLocalRank == 0):
            train_data_list.append(oPsm.feature_list)
            train_label_list.append(positive_int)
        elif oPsm.RealLabel == LabelRevTrain: # and (oPsm.iLocalRank == 0 or oPsm.iLocalRank == 1):
            train_data_list.append(oPsm.feature_list)
            train_label_list.append(negative_int)
                
    sys.stdout.write(str(len(train_data_list)) + "\t")
    
    train_data_np = np.array(train_data_list)[:, feature_selection_list]
    train_label_np = np.array(train_label_list)
    
    num_positive = float((train_label_np==LabelPositive).sum())
    num_negative = float((train_label_np==LabelNegative).sum())
    class_weight_dict = {0: (num_positive/(num_negative+num_positive)), 1:(num_negative/(num_negative+num_positive))}
   
    clf = Sequential()
    
    clf.add(Dense(32, input_dim=len(feature_selection_list), init='uniform', activation='sigmoid'))
    clf.add(Dropout(0))
    clf.add(Dense(32, init='uniform', activation='sigmoid'))
    clf.add(Dropout(0))
    clf.add(Dense(1, init='uniform', activation='sigmoid'))
    
    # sgd = SGD(lr=0.1, decay=1e-6, momentum=0.9, nesterov=True)
    clf.compile(loss='binary_crossentropy',
              optimizer='rmsprop',
              metrics=['accuracy'])
    # clf.compile(loss='categorical_crossentropy', optimizer=SGD(lr=0.01, momentum=0.9, nesterov=True))
    
    clf.fit(train_data_np, train_label_np, batch_size=32, nb_epoch=2, verbose=1)#,class_weight=class_weight_dict)
    #predict_np = logreg.predict(train_data_np)

    a = clf.get_weights()
    b=np.dot(a[0], np.dot(a[2], a[4]))
    sum_float = 0
    for val in b:
        sum_float += abs(val)
    for i in range(len(feature_selection_list)):
        sys.stdout.write('%.3f\n' % abs(b[i]/sum_float))
    return



    # # test
    test_unknown_np = np.array(test_data_list)[:, feature_selection_list]
    # test_unknown_np = preprocessing.normalize(test_unknown_np, norm='max',axis=0)
    predict_np = clf.predict_proba(test_unknown_np)
    #predict_np = logreg.predict_proba(train_data_np)
    for i in range(len(predict_np)):
        psm_rank_list[i].fPredictProbability = predict_np[i]
        psm_rank_list[i].ML_feature.append(predict_np[i])
        # psm_list[id_list[i]].ML_feature.append(predict_np[i])
    # fdr_rank(psm_list, 3) 
    # return None
        
    psm_new_list = re_rank(psm_rank_list)
    # np.savetxt("predict_probability.txt", predict_np)
        #print str(len(psm_list_selected))

    psm_filtered_list_local = show_Fdr(psm_new_list, None, fdr=fdr_given)
    
    psm_filtered_list = []
    # show_Fdr_charge(psm_list)
    psm_filtered_list.extend(psm_filtered_list_local)
    
    return psm_filtered_list
    
    a = clf.get_weights()
    
    # for i_int in range(len(a)):
        # np.savetxt("/media/xgo/Seagate/Proteomics/Experiments/Benchmark/Sipros10/D10/PsmLatest/Table2/deeplearning_weight_"+str(i_int)+".txt", a[i_int], delimiter = '\t')
    # np.savetxt("/media/xgo/Seagate/Proteomics/Experiments/Benchmark/Sipros10/D10/PsmLatest/Table2/deeplearning_weight_2.txt", a[1], delimiter = '\t')
    # np.savetxt("/media/xgo/Seagate/Proteomics/Experiments/Benchmark/Sipros10/D10/PsmLatest/Table2/deeplearning_weight_3.txt", a[2], delimiter = '\t')
    print(str(len(a)))
    b=np.dot(a[0], np.dot(a[2], a[4]))
    sum_float = 0
    for val in b:
        sum_float += abs(val)
    idx = 0
    for i in range(len(feature_name_list)):
        if i in feature_selection_list:
            sys.stdout.write('%.5f' % abs(b[idx]/sum_float))
            idx += 1
        sys.stdout.write('\t')
    order_list = [2, 0, 1, 6, 4, 5, 12, 3, 7, 8, 9, 10, 11]
    print("\nfeature importance:")
    for ind in order_list:
        sys.stdout.write('%.5f\n' % abs(b[ind]/sum_float))
    return psm_filtered_list


from sklearn.ensemble import AdaBoostClassifier

def test_AdaBoost(psm_dict, psm_list, psm_neg_list, fdr_given=None, psm_left_list=None):
    # machine learning
    # # construct training data
    psm_filtered_list = []
    #psm_list_selected = []
    train_data_list = []
    train_label_list = []
    test_data_list = []
    positive_int = 1
    negative_int = 0
    psm_rank_list = [] 
    for oPsm in psm_list:
        if len(oPsm.feature_list) == 0:
            print('check')
        if oPsm.TrainingLabel == LabelFiltered:
            continue
        if oPsm.RealLabel == LabelRevReserve:
            # train_data_list.append(oPsm.feature_list)
            # train_label_list.append(negative_int)
            continue
        test_data_list.append(oPsm.feature_list)
        psm_rank_list.append(oPsm)
        if oPsm.RealLabel != LabelRevTrain and (oPsm.iLocalRank == 0):
            train_data_list.append(oPsm.feature_list)
            train_label_list.append(positive_int)
        elif oPsm.RealLabel == LabelRevTrain: # and (oPsm.iLocalRank == 0 or oPsm.iLocalRank == 1):
            train_data_list.append(oPsm.feature_list)
            train_label_list.append(negative_int)
    
    sys.stdout.write(str(len(train_data_list)) + "\t")
    
    train_data_np = np.array(train_data_list)[:, feature_selection_list]
    train_label_np = np.array(train_label_list)
    
    clf = AdaBoostClassifier(# base_estimator=DecisionTreeClassifier(min_samples_split=800,min_samples_leaf=100),
                             n_estimators=200,
                             random_state=50)
    
    clf.fit(train_data_np, train_label_np)
    #predict_np = logreg.predict(train_data_np)

    for i in range(len(feature_selection_list)):
        sys.stdout.write('%.3f\n' % clf.feature_importances_[i])
    return


    # # test
    test_unknown_np = np.array(test_data_list)[:, feature_selection_list]
    predict_np = clf.predict_proba(test_unknown_np)
    #predict_np = logreg.predict_proba(train_data_np)
    for i in range(len(predict_np)):
        psm_rank_list[i].fPredictProbability = predict_np[i, 1]
        psm_rank_list[i].ML_feature.append(predict_np[i, 1])
    # fdr_rank(psm_list, 2)
    # return None
        
    psm_new_list = re_rank(psm_rank_list)
    # del psm_list[:]
    # psm_list.extend(psm_new_list)
    # np.savetxt("predict_probability.txt", predict_np)
        #print str(len(psm_list_selected))

    psm_filtered_list_local = show_Fdr(psm_new_list, None, fdr=fdr_given)
    # show_Fdr_charge(psm_list)
    psm_filtered_list.extend(psm_filtered_list_local)
    
    return psm_filtered_list

    idx = 0
    for i in range(len(feature_name_list)):
        if i in feature_selection_list:
            sys.stdout.write('%.5f' % clf.feature_importances_[idx])
            idx += 1
        sys.stdout.write('\t')
    
    sys.stdout.write('\n')
    
    order_list = [2, 0, 1, 6, 4, 5, 12, 3, 7, 8, 9, 10, 11]
    print("\nfeature importance:")
    for ind in order_list:
        sys.stdout.write('%.5f\n' % clf.feature_importances_[ind])
    
    return psm_filtered_list
    
    
from sklearn.svm import SVC

def test_SVM(psm_dict, psm_list, psm_neg_list, fdr_given=None, psm_left_list=None):
    # machine learning
    # # construct training data
    psm_filtered_list = []
    #psm_list_selected = []
    data_list = []
    label_list = []
    unknown_list = []
    id_list = [] 
    for oPsm in psm_list:
        if len(oPsm.feature_list) == 0:
            print('check')
        if oPsm.TrainingLabel == LabelFiltered:
            continue
        unknown_list.append(oPsm.feature_list)
        id_list.append(oPsm.iInnerId)
        if oPsm.TrainingLabel != LabelUnknown:
            data_list.append(oPsm.feature_list)
            label_list.append(oPsm.TrainingLabel)
    
           
    for oPsm in psm_neg_list:
        data_list.append(oPsm.feature_list)
        label_list.append(oPsm.TrainingLabel)
                
    data_np = np.array(data_list)
    sys.stdout.write(str(len(data_list)) + "\t")
    label_np = np.array(label_list)
    train_data_np = data_np[:, feature_selection_list]
    
    clf = SVC(random_state=50,
              kernel='linear',
              class_weight='balanced',
              cache_size=4000,
              C=1,
              gamma='auto')
    clf.fit(train_data_np, label_np)
    #predict_np = logreg.predict(train_data_np)

    # # test
    unknown_np = np.array(unknown_list)
    test_unknown_np = unknown_np[:, feature_selection_list]
    predict_np = clf.decision_function(test_unknown_np)
    #predict_np = logreg.predict_proba(train_data_np)
    for i in range(len(predict_np)):
        psm_list[id_list[i]].fPredictProbability = predict_np[i]
        
    psm_new_list = re_rank(psm_list)
    del psm_list[:]
    psm_list.extend(psm_new_list)
    # np.savetxt("predict_probability.txt", predict_np)
        #print str(len(psm_list_selected))
    if not psm_left_list is None:
        psm_filtered_list_local = show_Fdr(psm_left_list, None, fdr=fdr_given)
    else:
        psm_filtered_list_local = show_Fdr(psm_list, None, fdr=fdr_given)
    # show_Fdr_charge(psm_list)
    psm_filtered_list.extend(psm_filtered_list_local)
    '''
    idx = 0
    for i in range(len(feature_name_list)):
        if i in feature_selection_list:
            sys.stdout.write('%.3f' % clf.feature_importances_[idx])
            idx += 1
        sys.stdout.write('\t')
    '''    
    sys.stdout.write('\n')
        
    return psm_filtered_list

from sklearn.naive_bayes import GaussianNB

def test_naive_bayes(psm_dict, psm_list, psm_neg_list, fdr_given=None, psm_left_list=None):
    # machine learning
    # # construct training data
    psm_filtered_list = []
    #psm_list_selected = []
    data_list = []
    label_list = []
    unknown_list = []
    id_list = [] 
    for oPsm in psm_list:
        if len(oPsm.feature_list) == 0:
            print('check')
        if oPsm.TrainingLabel == LabelFiltered:
            continue
        unknown_list.append(oPsm.feature_list)
        id_list.append(oPsm.iInnerId)
        if oPsm.TrainingLabel != LabelUnknown:
            data_list.append(oPsm.feature_list)
            label_list.append(oPsm.TrainingLabel)
    
           
    for oPsm in psm_neg_list:
        data_list.append(oPsm.feature_list)
        label_list.append(oPsm.TrainingLabel)
                
    data_np = np.array(data_list)
    sys.stdout.write(str(len(data_list)) + "\t")
    label_np = np.array(label_list)
    train_data_np = data_np[:, feature_selection_list]
    train_data_np = preprocessing.normalize(train_data_np, norm='l2',axis=0)
    
    clf = GaussianNB()
    clf.fit(train_data_np, label_np)
    #predict_np = logreg.predict(train_data_np)

    # # test
    unknown_np = np.array(unknown_list)
    test_unknown_np = unknown_np[:, feature_selection_list]
    test_unknown_np = preprocessing.normalize(test_unknown_np, norm='l2',axis=0)
    
    predict_np = clf.predict_proba(test_unknown_np)
    #predict_np = logreg.predict_proba(train_data_np)
    for i in range(len(predict_np)):
        psm_list[id_list[i]].fPredictProbability = predict_np[i, 1]
        
    psm_new_list = re_rank(psm_list)
    del psm_list[:]
    psm_list.extend(psm_new_list)
    # np.savetxt("predict_probability.txt", predict_np)
        #print str(len(psm_list_selected))
    if not psm_left_list is None:
        psm_filtered_list_local = show_Fdr(psm_left_list, None, fdr=fdr_given)
    else:
        psm_filtered_list_local = show_Fdr(psm_list, None, fdr=fdr_given)
    # show_Fdr_charge(psm_list)
    psm_filtered_list.extend(psm_filtered_list_local)
    '''
    idx = 0
    for i in range(len(feature_name_list)):
        if i in feature_selection_list:
            sys.stdout.write('%.3f' % clf.feature_importances_[idx])
            idx += 1
        sys.stdout.write('\t')
    '''    
    sys.stdout.write('\n')
        
    return psm_filtered_list

def test_random_forest(psm_dict, psm_list, psm_neg_list, fdr_given=None, psm_left_list=None, output_folder=None, input_file=None):
    # machine learning
    train_data_list = []
    train_label_list = []
    test_data_list = []
    num_feature_int = (31 + len(ptm_selection_list))
    positive_int = 1
    negative_int = 0
    psm_rank_list = []         
    for oPsm in psm_list:
        if len(oPsm.feature_list) != num_feature_int:
            pass
            # print 'check'
        if oPsm.TrainingLabel == LabelFiltered:
            continue
        if oPsm.RealLabel == LabelRevReserve:
            # train_data_list.append(oPsm.feature_list)
            # train_label_list.append(negative_int)
            continue
        # if oPsm.iLocalRank == 0 or oPsm.iLocalRank == 1:
        test_data_list.append(oPsm.feature_list)
        psm_rank_list.append(oPsm)
        # if oPsm.TrainingLabel != LabelUnknown and oPsm.iLocalRank == 0:
        if oPsm.RealLabel != LabelRevTrain and (oPsm.iLocalRank == 0):
            train_data_list.append(oPsm.feature_list)
            train_label_list.append(positive_int)
        elif oPsm.RealLabel == LabelRevTrain: # and (oPsm.iLocalRank == 0 or oPsm.iLocalRank == 1):
            train_data_list.append(oPsm.feature_list)
            train_label_list.append(negative_int)
                
    sys.stdout.write(str(len(train_data_list)) + "\t")
        
    train_data_np = np.array(train_data_list)[:, feature_selection_list]
    train_label_np = np.array(train_label_list)
    
     # only forward left
    unique_np = np.unique(train_label_np)
    if unique_np.shape[0] == 1:
        psm_filtered_list_local = show_Fdr(psm_list, None, fdr=fdr_given)
        return psm_filtered_list_local
    # train_data_np = preprocessing.normalize(train_data_np, norm='max',axis=0)
    
    clf = RandomForestClassifier(n_estimators=200, 
                                 n_jobs=-1, 
                                 class_weight='balanced', 
                                 oob_score=True,
                                 bootstrap=True,
                                 random_state=50,
                                 criterion="gini",
                                 max_features="auto",
                                 min_samples_leaf=50,
                                 min_samples_split=800)
    
    clf.fit(train_data_np, train_label_np)
    
    for i in range(len(feature_selection_list)):
        sys.stdout.write('%.5f' % clf.feature_importances_[i])
        sys.stdout.write('\n')
    return
    
    # # test
    test_unknown_np = np.array(test_data_list)[:, feature_selection_list]
    predict_np = clf.predict_proba(test_unknown_np)
    
    for i in range(len(predict_np)):
        psm_rank_list[i].fPredictProbability = predict_np[i, 1]
        psm_rank_list[i].ML_feature.append(predict_np[i, 1])
    # fdr_rank(psm_list, 1)
    # return None
    
    psm_new_list = re_rank(psm_rank_list)
    # del psm_list[:]
    # psm_list.extend(psm_new_list)
    '''
    for fdr_f in [0.001, 0.0025, 0.005, 0.0075, 0.01, 0.02]:
        temp_list = show_Fdr_varied(psm_list, fdr_f/3.0)
        folder_str = output_folder + 'fdr_' + str(fdr_f) +'/'
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
        generate_psm_pep_txt(input_file, folder_str, temp_list)
    '''

    psm_filtered_list_local = show_Fdr(psm_new_list, None, fdr=fdr_given)
    
    psm_filtered_list = []
    # show_Fdr_charge(psm_list)
    psm_filtered_list.extend(psm_filtered_list_local)
    
    return psm_filtered_list

    idx = 0
    for i in range(len(feature_name_list)):
        if i in feature_selection_list:
            sys.stdout.write('%.5f' % clf.feature_importances_[idx])
            idx += 1
        sys.stdout.write('\t')
    sys.stdout.write('\n')
    '''
    order_list = [2, 0, 1, 6, 4, 5, 12, 3, 7, 8, 9, 10, 11]
    print "\nfeature importance:"
    for ind in order_list:
        sys.stdout.write('%.5f\n' % clf.feature_importances_[ind])
    '''
    return psm_filtered_list

def test_stacking_MOE(psm_list):
    # # construct training data
    psm_filtered_list = []
    data_list = []
    label_list = []
    id_list = []
    D_list = []
    U_list = []
    D_id_list = []
    U_id_list = []     
    for oPsm in psm_list:
        if len(oPsm.feature_list) == 0:
            print('check')
        id_list.append(oPsm.iInnerId)
        if oPsm.TrainingLabel != LabelUnknown:
            data_list.append(oPsm.feature_list)
            label_list.append(oPsm.TrainingLabel)
        if oPsm.TrainingLabel == LabelNegative:
            D_list.append(oPsm.feature_list)
            D_id_list.append(oPsm.iInnerId)
        else:
            U_list.append(oPsm.feature_list)
            U_id_list.append(oPsm.iInnerId)
                
    data_np = np.array(data_list)
    sys.stdout.write(str(len(data_list)) + "\n")
    label_np = np.array(label_list)
    
    feature_list_local = [1, 2, 3, 16]
    D_np = np.array(D_list, dtype=np.float64)
    U_np = np.array(U_list, dtype=np.float64)
    
    # features and models
    feature_matrix = [[1, 2, 3, 15, 16, 17],    # score section
                      [5, 29],                  # mass section
                      [24],                     # digestion section
                      [30],                     # PTM section
                      [25, 26, 27, 28]]         # pep pro support section
    tier_1 = []
    D_phix = []
    U_phix = []
    for i in range(len(feature_matrix)):
        tier_1.append(linear_model.LogisticRegression(penalty='l2', 
                                             dual=False, 
                                             C=1, 
                                             fit_intercept=True, 
                                             intercept_scaling=1, 
                                             class_weight='balanced', 
                                             random_state=None, 
                                             solver='liblinear', 
                                             max_iter=100, 
                                             multi_class='ovr', 
                                             verbose=0, 
                                             warm_start=False, 
                                             n_jobs=-1))
        tier_1[i].fit(data_np[:, feature_matrix[i]], label_np)
        predict_np = tier_1[i].predict_proba(D_np[:, feature_matrix[i]])
        D_phix.append(predict_np[:, 0])
        predict_np = tier_1[i].predict_proba(U_np[:, feature_matrix[i]])
        U_phix.append(predict_np[:, 1])
    D_phi_np = np.transpose(np.array(D_phix))
    U_phi_np = np.transpose(np.array(U_phix))
    D_np = D_np[:, feature_list_local]
    U_np = U_np[:, feature_list_local]
    (D_pyx, U_pyx) = test_MoE_semi(D_np, U_np, D_phi_np, U_phi_np)
    # D_pyx = D_phi_np[:,1]
    # U_pyx = U_phi_np[:,1]
    
    for i in range(len(D_id_list)):
        psm_list[D_id_list[i]].fPredictProbability = 1 - D_pyx[i]
    
    for i in range(len(U_id_list)):
        psm_list[U_id_list[i]].fPredictProbability = U_pyx[i]
    
    psm_new_list = re_rank(psm_list)
    del psm_list[:]
    psm_list.extend(psm_new_list)
    psm_filtered_list_local = show_Fdr(psm_list, None, fdr=0.01)
    psm_filtered_list.extend(psm_filtered_list_local)
    sys.stdout.write('\n')
    
    return psm_filtered_list

def test_MoE_super(D_np, F_np, D_phi_np, F_phi_np):
    # normalize 
    preprocessing.normalize(D_np, axis=0,copy=False, norm='max')
    preprocessing.normalize(F_np, axis=0,copy=False, norm='max')
    # preprocessing.normalize(D_phi_np, axis=0,copy=False, norm='max')
    # preprocessing.normalize(F_phi_np, axis=0,copy=False, norm='max')
    # n sample #, m expert #, d feature #
    # initialize parameters: alpha, mu, sigma_square
    numD = np.float128(D_np.shape[0])
    numF = np.float128(F_np.shape[0])
    num = numD + numF
    d = np.float128(D_np.shape[1])
    m = np.float128(D_phi_np.shape[1])
    alpha = np.zeros(m, dtype=np.float128) + 1.0/float(m)
    mu_D = np.average(D_np, axis=0)
    mu_U = np.average(F_np, axis=0)
    mu_temp = (mu_D + mu_U)/2.0
    mu = np.tile(mu_temp, (m, 1))
    sigma_square = np.zeros(m, dtype=np.float128) + 1.0
    pi = math.pi
    
    # alpha = p(z,y|theta) 
    alpha_old = np.copy(alpha)
    mu_old = np.copy(mu)
    sigma_square_old = np.copy(sigma_square)
    
    # the iteration of E-M algorithm
    for _ite in range(100):
        # E-step
        D_minus_mu = []
        F_minus_mu = []
        for mu_j in mu:
            D_minus_mu.append(np.sum(np.multiply((D_np - mu_j), (D_np - mu_j)), axis=1)) # (1, n)
            F_minus_mu.append(np.sum(np.multiply((F_np - mu_j), (F_np - mu_j)), axis=1))
        D_minus_mu = np.array(D_minus_mu, dtype=np.float128) # (m, n)
        D_minus_mu = np.transpose(D_minus_mu) # (n, m)
        F_minus_mu = np.array(F_minus_mu, dtype=np.float128) # (m, n)
        F_minus_mu = np.transpose(F_minus_mu) # (n, m)
        
        D_px = np.multiply((np.power(2*pi*sigma_square, -d/2)), np.exp(np.divide(D_minus_mu, (-2*sigma_square))))
        D_alpha_px = np.multiply(alpha, D_px) # (n, m)
        D_phi_alpha_px = np.multiply(D_alpha_px , D_phi_np) # (n, m)
        D_phi_alpha_px_sum = np.sum(D_phi_alpha_px, axis=1)
        # # Q(z=j)
        D_hzj = np.divide(D_phi_alpha_px , D_phi_alpha_px_sum[:, None]) # (n, m)
        
        F_px = np.multiply((np.power(2*pi*sigma_square, -d/2)), np.exp(np.divide(F_minus_mu, (-2*sigma_square))))
        F_alpha_px = np.multiply(alpha, F_px) # (n, m)
        F_phi_alpha_px = np.multiply(F_alpha_px , F_phi_np) # (n, m)
        F_phi_alpha_px_sum = np.sum(F_phi_alpha_px, axis=1)
        # # Q(z=j)
        F_hzj = F_phi_alpha_px / F_phi_alpha_px_sum[:,None] # (n, m)
        
        # log-likelihood
        part1 = np.sum(np.log(D_phi_alpha_px_sum), axis=0)
        part2 = np.sum(np.log(F_phi_alpha_px_sum), axis=0)
        ftheta = part1 + part2
        print(str(ftheta))
        
        # M -step
        # # alpha
        alpha = np.divide((np.sum(D_hzj, axis=0) + np.sum(F_hzj, axis=0)), (num))
        # # mu_j, sigma_j
        for j in range(int(m)):
            # sigma square
            sigma_square[j] = (np.sum(D_hzj[:,j]*np.sum(np.multiply((D_np - mu[j]), (D_np - mu[j])), axis = 1)) + 
                               np.sum(F_hzj[:,j]*np.sum(np.multiply((F_np - mu[j]), (F_np - mu[j])), axis = 1)))
            sigma_square[j] /= (d*alpha[j] * num)
            # mu
            mu[j] = (np.sum(D_np*((D_hzj[:,j])[:,None]), axis=0)+np.sum(F_np*((F_hzj[:,j])[:,None]), axis=0))/(alpha[j] * num)
        
        # check if converge
        error_sum = np.sum(np.absolute(mu-mu_old)) + np.sum(np.absolute(alpha-alpha_old)) + np.sum(np.absolute(sigma_square-sigma_square_old))
        error_sum /= float(mu.size + alpha.size + sigma_square.size)
        alpha_old = np.copy(alpha)
        mu_old = np.copy(mu)
        sigma_square_old = np.copy(sigma_square)
        
        if error_sum < 0.0001:
            print('Converged:\t%d,\tError:\t%f' % (_ite, error_sum))
            '''
            print alpha
            print mu
            print sigma_square
            '''
            break
        else:
            print('Iteration:\t%d,\tError:\t%f' % (_ite, error_sum))
            '''
            print alpha
            print mu
            print sigma_square
            '''
    
    # assign probability
    D_minus_mu = []
    F_minus_mu = []
    for mu_j in mu:
        D_minus_mu.append(np.sum((D_np - mu_j)*(D_np - mu_j), axis=1)) # (1, n)
        F_minus_mu.append(np.sum((F_np - mu_j)*(F_np - mu_j), axis=1))
    D_minus_mu = np.array(D_minus_mu, dtype=np.float128) # (m, n)
    D_minus_mu = np.transpose(D_minus_mu) # (n, m)
    F_minus_mu = np.array(F_minus_mu, dtype=np.float128) # (m, n)
    F_minus_mu = np.transpose(F_minus_mu) # (n, m)
    D_px = np.multiply((np.power(2*pi*sigma_square, -d/2)), np.exp(np.divide(D_minus_mu, (-2*sigma_square))))
    D_alpha_px = np.multiply(alpha, D_px) # (n, m)
    D_alpha_px_sum = np.sum(D_alpha_px, axis=1)
    D_gx = D_alpha_px / D_alpha_px_sum[:,None]
    D_phi_gx = D_phi_np * D_gx
    D_pyx = np.sum(D_phi_gx, axis=1)
    
    F_px = np.multiply((np.power(2*pi*sigma_square, -d/2)), np.exp(np.divide(F_minus_mu, (-2*sigma_square))))
    F_alpha_px = np.multiply(alpha, F_px) # (n, m)
    F_alpha_px_sum = np.sum(F_alpha_px, axis=1)
    F_gx = F_alpha_px / F_alpha_px_sum[:, None]
    F_phi_gx = F_phi_np * F_gx
    U_pyx = np.sum(F_phi_gx, axis=1)
    
    return (D_pyx, U_pyx)

def test_MoE_semi(D_np, F_np, D_phi_np, F_phi_np):
    # normalize 
    preprocessing.normalize(D_np, axis=0,copy=False, norm='max')
    preprocessing.normalize(F_np, axis=0,copy=False, norm='max')
    F_phi_2_np = F_phi_np
    F_phi_1_np = 1.0 - F_phi_np 
    # preprocessing.normalize(D_phi_np, axis=0,copy=False, norm='max')
    # preprocessing.normalize(F_phi_np, axis=0,copy=False, norm='max')
    # n sample #, m expert #, d feature #
    # initialize parameters: alpha_j1, mu_j1, sigma_square_j1
    numD = np.float128(D_np.shape[0])
    numF = np.float128(F_np.shape[0])
    num = numD + numF
    d = np.float128(D_np.shape[1])
    m = np.float128(D_phi_np.shape[1])
    alpha_j1 = np.zeros(m, dtype=np.float128) + 1.0/float(m)
    alpha_j2 = np.zeros(m, dtype=np.float128) + 1.0/float(m)
    mu_D = np.average(D_np, axis=0)
    mu_U = np.average(F_np, axis=0)
    mu_temp = (mu_D + mu_U)/2.0
    mu_j1 = np.tile(mu_temp, (m, 1))
    mu_j2 = np.tile(mu_temp, (m, 1))
    sigma_square_j1 = np.zeros(m, dtype=np.float128) + 1.0
    sigma_square_j2 = np.zeros(m, dtype=np.float128) + 1.0
    pi = math.pi
    
    # alpha = p(z,y|theta) 
    alpha_j1_old = np.copy(alpha_j1)
    alpha_j2_old = np.copy(alpha_j2)
    mu_j1_old = np.copy(mu_j1)
    mu_j2_old = np.copy(mu_j2)
    sigma_square_j1_old = np.copy(sigma_square_j1)
    sigma_square_j2_old = np.copy(sigma_square_j2)
    
    # the iteration of E-M algorithm
    for _ite in range(100):
        # E-step
        D_minus_mu_j1 = []
        F_minus_mu_j1 = []
        F_minus_mu_j2 = []
        for mu_j in mu_j1:
            D_minus_mu_j1.append(np.sum(np.multiply((D_np - mu_j), (D_np - mu_j)), axis=1)) # (1, n)
            F_minus_mu_j1.append(np.sum(np.multiply((F_np - mu_j), (F_np - mu_j)), axis=1))
        for mu_j in mu_j2:
            F_minus_mu_j2.append(np.sum(np.multiply((F_np - mu_j), (F_np - mu_j)), axis=1))
        D_minus_mu_j1 = np.array(D_minus_mu_j1, dtype=np.float128) # (m, n)
        D_minus_mu_j1 = np.transpose(D_minus_mu_j1) # (n, m)
        F_minus_mu_j1 = np.array(F_minus_mu_j1, dtype=np.float128) # (m, n)
        F_minus_mu_j1 = np.transpose(F_minus_mu_j1) # (n, m)
        F_minus_mu_j2 = np.array(F_minus_mu_j2, dtype=np.float128) # (m, n)
        F_minus_mu_j2 = np.transpose(F_minus_mu_j2) # (n, m)
        D_px_zy1 = np.multiply((np.power(2*pi*sigma_square_j1, -d/2)), np.exp(np.divide(D_minus_mu_j1, (-2*sigma_square_j1))))
        D_alpha_px_zy1 = np.multiply(alpha_j1, D_px_zy1) # (n, m)
        D_phi_alpha_px_zy1 = np.multiply(D_alpha_px_zy1 , D_phi_np) # (n, m)
        D_phi_alpha_px_zy1_sum = np.sum(D_phi_alpha_px_zy1, axis=1)
        # # Q(z=j)
        Qzj = np.divide(D_phi_alpha_px_zy1 , D_phi_alpha_px_zy1_sum[:, None]) # (n, m)
        # D_alphaPx_log = np.exp(np.log(alpha_j1) + np.log(pow(2*pi*sigma_square_j1, d/2)) + (D_minus_mu_j1*(-1)/(2*sigma_square_j1))) # (n, m)
        F_px_zy1 = np.multiply((np.power(2*pi*sigma_square_j1, -d/2)), np.exp(np.divide(F_minus_mu_j1, (-2*sigma_square_j1))))
        F_px_zy2 = np.multiply((np.power(2*pi*sigma_square_j2, -d/2)), np.exp(np.divide(F_minus_mu_j2, (-2*sigma_square_j2))))
        F_alpha_px_zy1 = np.multiply(alpha_j1, F_px_zy1) # (n, m)
        F_alpha_px_zy2 = np.multiply(alpha_j2, F_px_zy2) # (n, m)
        F_phi_alpha_px_zy1 = np.multiply(F_alpha_px_zy1 , F_phi_1_np) # (n, m)
        F_phi_alpha_px_zy2 = np.multiply(F_alpha_px_zy2 , F_phi_2_np) # (n, m)
        F_phi_alpha_px_zy1_sum = np.sum(F_phi_alpha_px_zy1, axis=1)
        F_phi_alpha_px_zy2_sum = np.sum(F_phi_alpha_px_zy2, axis=1)
        F_phi_alpha_px_zy_sum = F_phi_alpha_px_zy1_sum + F_phi_alpha_px_zy2_sum
        Qzjy1 = F_phi_alpha_px_zy1 / F_phi_alpha_px_zy_sum[:,None] # (n, m)
        Qzjy2 = F_phi_alpha_px_zy2 / F_phi_alpha_px_zy_sum[:,None] # (n, m)
        
        # log-likelihood
        part1 = np.sum(np.log(D_phi_alpha_px_zy1_sum), axis=0)
        part2 = np.sum(np.log(F_phi_alpha_px_zy_sum), axis=0)
        ftheta = part1 + part2
        print(str(ftheta))
        
        # M -step
        # # alpha
        alpha_j1 = np.divide((np.sum(Qzj, axis=0) + np.sum(Qzjy1, axis=0)), (num))
        alpha_j2 = np.divide(np.sum(Qzjy2, axis=0), (num))
        # # mu_j, sigma_j
        for j in range(int(m)):
            # sigma square
            sigma_square_j1[j] = (np.sum(Qzj[:,j]*np.sum(np.multiply((D_np - mu_j1[j]), (D_np - mu_j1[j])), axis=1)) + 
                               np.sum(Qzjy1[:,j]*np.sum(np.multiply((F_np - mu_j1[j]), (F_np - mu_j1[j])), axis = 1)))
            sigma_square_j1[j] /= (d*alpha_j1[j] * num)
            sigma_square_j2[j] = np.sum(Qzjy2[:,j]*np.sum(np.multiply((F_np - mu_j2[j]), (F_np - mu_j2[j])), axis = 1))
            sigma_square_j2[j] /= (d*alpha_j2[j] * num)
            # mu
            mu_j1[j] = (np.sum(D_np*((Qzj[:,j])[:,None]), axis=0)+np.sum(F_np*((Qzjy1[:,j])[:,None]), axis=0))/(alpha_j1[j] * num)
            mu_j2[j] = (np.sum(F_np*((Qzjy2[:,j])[:,None]), axis=0))/(alpha_j2[j] * num)
        
        # check if converge
        error_sum = np.sum(np.absolute(mu_j1-mu_j1_old)) + np.sum(np.absolute(alpha_j1-alpha_j1_old)) + np.sum(np.absolute(sigma_square_j1-sigma_square_j1_old))
        error_sum += np.sum(np.absolute(mu_j2-mu_j2_old)) + np.sum(np.absolute(alpha_j2-alpha_j2_old)) + np.sum(np.absolute(sigma_square_j2-sigma_square_j2_old))
        error_sum /= 2.0*float(mu_j1.size + alpha_j1.size + sigma_square_j1.size)
        alpha_j1_old = np.copy(alpha_j1)
        mu_j1_old = np.copy(mu_j1)
        sigma_square_j1_old = np.copy(sigma_square_j1)
        alpha_j2_old = np.copy(alpha_j2)
        mu_j2_old = np.copy(mu_j2)
        sigma_square_j2_old = np.copy(sigma_square_j2)
        
        if error_sum < 0.0001:
            print('Converged:\t%d,\tError:\t%f' % (_ite, error_sum))
            '''
            print alpha_j1
            print mu_j1
            print sigma_square_j1
            '''
            break
        else:
            print('Iteration:\t%d,\tError:\t%f' % (_ite, error_sum))
            '''
            print alpha_j1
            print mu_j1
            print sigma_square_j1
            '''
    
    # assign probability
    D_minus_mu_j1 = []
    F_minus_mu_j1 = []
    F_minus_mu_j2 = []
    for mu_j in mu_j1:
        D_minus_mu_j1.append(np.sum((D_np - mu_j)*(D_np - mu_j), axis=1)) # (1, n)
        F_minus_mu_j1.append(np.sum((F_np - mu_j)*(F_np - mu_j), axis=1))
    for mu_j in mu_j2:
        F_minus_mu_j2.append(np.sum(np.multiply((F_np - mu_j), (F_np - mu_j)), axis=1))
    D_minus_mu_j1 = np.array(D_minus_mu_j1, dtype=np.float128) # (m, n)
    D_minus_mu_j1 = np.transpose(D_minus_mu_j1) # (n, m)
    F_minus_mu_j1 = np.array(F_minus_mu_j1, dtype=np.float128) # (m, n)
    F_minus_mu_j1 = np.transpose(F_minus_mu_j1) # (n, m)
    F_minus_mu_j2 = np.array(F_minus_mu_j2, dtype=np.float128) # (m, n)
    F_minus_mu_j2 = np.transpose(F_minus_mu_j2) # (n, m)
    D_px_zy1 = np.multiply((np.power(2*pi*sigma_square_j1, -d/2)), np.exp(np.divide(D_minus_mu_j1, (-2*sigma_square_j1))))
    D_alpha_px_zy1 = np.multiply(alpha_j1, D_px_zy1) # (n, m)
    D_alpha_px_zy1_sum = np.sum(D_alpha_px_zy1, axis=1)
    D_pzj_x = D_alpha_px_zy1 / D_alpha_px_zy1_sum[:,None]
    D_phi_pzj_x = D_phi_np * D_pzj_x
    D_pyx = np.sum(D_phi_pzj_x, axis=1)
    
    F_px_zy1 = np.multiply((np.power(2*pi*sigma_square_j1, -d/2)), np.exp(np.divide(F_minus_mu_j1, (-2*sigma_square_j1))))
    F_px_zy2 = np.multiply((np.power(2*pi*sigma_square_j2, -d/2)), np.exp(np.divide(F_minus_mu_j2, (-2*sigma_square_j2))))
    F_alpha_px_zy1 = np.multiply(alpha_j1, F_px_zy1) # (n, m)
    F_alpha_px_zy2 = np.multiply(alpha_j2, F_px_zy2) # (n, m)
    F_alpha_px_zy = F_alpha_px_zy1 + F_alpha_px_zy2
    F_alpha_px_zy_sum = np.sum(F_alpha_px_zy, axis=1)
    F_pzj_x = F_alpha_px_zy / F_alpha_px_zy_sum[:, None]
    F_phi_pzj_x = F_phi_2_np * F_pzj_x
    U_pyx = np.sum(F_phi_pzj_x, axis=1)
    
    return (D_pyx, U_pyx)

    
    
def logistic_regression(psm_dict, psm_list):
    # machine learning
    # # construct training data
    psm_filtered_list = []
    psm_list_selected = []
    for key, lPsm in psm_dict.items():
        sys.stdout.write(key + "\t")
        psm_list_selected = lPsm
        data_list = []
        label_list = []
        unknown_list = []
        id_list = []
        for oPsm in lPsm:
            #feature_list = []
            #oPsm.set_feature(feature_list)
            if len(oPsm.feature_list) == 0:
                print('check')
            unknown_list.append(oPsm.feature_list)
            id_list.append(oPsm.iInnerId)
            if oPsm.TrainingLabel != LabelUnknown:
                data_list.append(oPsm.feature_list)
                label_list.append(oPsm.TrainingLabel)

        data_np = np.array(data_list)
        sys.stdout.write(str(len(data_list)) + "\t")
        label_np = np.array(label_list)
        
        '''
        for i in range(len(feature_name_list)):
            del feature_selection_list[:]
            feature_selection_list.append(i)
            for j in range(len(feature_name_list)):
                if j != i:
                    feature_selection_list.append(j)
        '''    
        '''
        del feature_selection_list[:]
        for i in range(len(feature_name_list)):
            if i < 6 or i > 23 or ( i>=15 and i<=17):
                feature_selection_list.append(i)
        #feature_selection_list.extend([1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28])
        '''
        train_data_np = data_np[:, feature_selection_list]
    
    # # training
        class_dict = {0: 1, 1:1000}
        logreg = linear_model.LogisticRegression(penalty='l2', 
                                             dual=False, 
                                             C=1, 
                                             fit_intercept=True, 
                                             intercept_scaling=1, 
                                             class_weight='balanced', 
                                             random_state=None, 
                                             solver='liblinear', 
                                             max_iter=100, 
                                             multi_class='ovr', 
                                             verbose=0, 
                                             warm_start=False, 
                                             n_jobs=-1)
        logreg.fit(train_data_np, label_np)
        predict_np = logreg.predict(train_data_np)
        #np.savetxt("xxx.txt", train_data_np)
        #return

    #show_TP_TN_FP_FN(label_np, predict_np)
    #print 'Actual number of iterations for all classes.'
        #print logreg.n_iter_
    # # test
        unknown_np = np.array(unknown_list)
        test_unknown_np = unknown_np[:, feature_selection_list]
        predict_np = logreg.predict_proba(test_unknown_np)

        for i in range(len(predict_np)):
            psm_list[id_list[i]].fPredictProbability = predict_np[i, 1]

    # np.savetxt("predict_probability.txt", predict_np)
        #print str(len(psm_list_selected))
        psm_filtered_list_local = show_Fdr(psm_list_selected, key)
        psm_filtered_list.extend(psm_filtered_list_local)
    
    #print 'Coefficient of the features in the decision function:'
    #print logreg.coef_
        idx = 0
        for i in range(len(feature_name_list)):
            if i in feature_selection_list:
                sys.stdout.write('%.3f' % logreg.coef_[0][idx])
                idx += 1
            sys.stdout.write('\t')
        sys.stdout.write('\n')
        #print str(logreg.intercept_)
    '''
    for x in logreg.coef_[0]:
        sys.stdout.write(str(x))
        sys.stdout.write('\t')
    '''
    return psm_filtered_list

def fdr_rank(psm_list, feature_id):
    list_sorted = sorted(psm_list, key=lambda x: (x.ML_feature[feature_id]) , reverse=True)
    num_Fwd = 0
    num_Rev = 0
    num_Shu = 0
    for oPsm in list_sorted:
        if oPsm.RealLabel == LabelFwd:
            num_Fwd += 1
        elif oPsm.RealLabel == LabelNegative:
            num_Rev += 1
        else:
            num_Shu += 1
        (_FDR_accept, FDR_value) = FDR_calculator(num_Shu, num_Fwd)
        oPsm.FDR_feature.append(FDR_value)

def re_rank(psm_list):
    psm_new_list = []
    psm_dict = {}
    for oPsm in psm_list:
        sId = oPsm.FileName + '_' + str(oPsm.ScanNumber)
        if sId in psm_dict:
            if oPsm.fPredictProbability > psm_dict[sId].fPredictProbability:
                psm_dict[sId] = oPsm
        else:
            psm_dict[sId] = oPsm
    
    for _key, value in psm_dict.items():
        psm_new_list.append(value)
    
    return psm_new_list

def test_MoE_ensamble(psm_list):
    psm_filtered_list = []
    id_list = []
    D_list = []
    U_list = []
    D_id_list = []
    U_id_list = []
    D_phi_list = []
    U_phi_list = []     
    for oPsm in psm_list:
        id_list.append(oPsm.iInnerId)
        if oPsm.TrainingLabel == LabelNegative:
            D_list.append(oPsm.feature_list)
            D_id_list.append(oPsm.iInnerId)
            D_phi_list.append(oPsm.ML_feature)
        else:
            U_list.append(oPsm.feature_list)
            U_id_list.append(oPsm.iInnerId)
            U_phi_list.append(oPsm.ML_feature)

    D_np = np.array(D_list, dtype=np.float64)
    U_np = np.array(U_list, dtype=np.float64)
    
    
    D_phi_np = np.array(D_phi_list, dtype=np.float128)
    U_phi_np = np.array(U_phi_list, dtype=np.float128)
    feature_list_local = [1, 2, 3, 5, 15, 16, 17, 24, 25, 26, 27, 28, 29]
    D_np = D_np[:, feature_list_local]
    U_np = U_np[:, feature_list_local]
    # (D_pyx, U_pyx) = test_MoE_semi(D_np, U_np, D_phi_np, U_phi_np)
    (D_pyx, U_pyx) = test_MoE_super(D_np, U_np, D_phi_np, U_phi_np)
    # D_pyx = D_phi_np[:,1]
    # U_pyx = U_phi_np[:,1]
    
    for i in range(len(D_id_list)):
        psm_list[D_id_list[i]].fPredictProbability = 1 - D_pyx[i]
    
    for i in range(len(U_id_list)):
        psm_list[U_id_list[i]].fPredictProbability = U_pyx[i]
    
    psm_new_list = re_rank(psm_list)
    del psm_list[:]
    psm_list.extend(psm_new_list)
    psm_filtered_list_local = show_Fdr(psm_list, None, fdr=0.01)
    psm_filtered_list.extend(psm_filtered_list_local)
    sys.stdout.write('\n')
    
    return psm_filtered_list

def test_LR_ensamble(psm_dict, psm_list, psm_neg_list, fdr_given=None, psm_left_list=None, output_folder=None, input_file=None):
    # machine learning
    # # construct training data
    psm_filtered_list = []
    #psm_list_selected = []
    train_data_list = []
    train_label_list = []
    test_data_list = []   
    positive_int = 1
    negative_int = 0
    psm_rank_list = [] 
    for oPsm in psm_list:
        if oPsm.TrainingLabel == LabelFiltered:
            continue
        if oPsm.RealLabel == LabelRevReserve:
            # train_data_list.append(oPsm.feature_list)
            # train_label_list.append(negative_int)
            continue
        test_data_list.append(oPsm.ML_feature)
        psm_rank_list.append(oPsm)
        if oPsm.RealLabel != LabelRevTrain and (oPsm.iLocalRank == 0):
            train_data_list.append(oPsm.ML_feature)
            train_label_list.append(positive_int)
        elif oPsm.RealLabel == LabelRevTrain: # and (oPsm.iLocalRank == 0 or oPsm.iLocalRank == 1):
            train_data_list.append(oPsm.ML_feature)
            train_label_list.append(negative_int)
                
    sys.stdout.write(str(len(train_data_list)) + "\t")
    
    train_data_np = np.array(train_data_list)
    train_label_np = np.array(train_label_list)
    
    logreg = linear_model.LogisticRegression(penalty='l2', 
                                             dual=False, 
                                             C=1, 
                                             fit_intercept=True, 
                                             intercept_scaling=1, 
                                             # class_weight='balanced',
                                             # class_weight=class_weight_dict, 
                                             random_state=None, 
                                             solver='liblinear', 
                                             max_iter=100, 
                                             multi_class='ovr', 
                                             verbose=0, 
                                             warm_start=False, 
                                             n_jobs=-1)
    logreg.fit(train_data_np, train_label_np)
    #predict_np = logreg.predict(train_data_np)

    # # test
    test_unknown_np = np.array(test_data_list)
    predict_np = logreg.predict_proba(test_unknown_np)
    #predict_np = logreg.predict_proba(train_data_np)
    for i in range(len(predict_np)):
        psm_rank_list[i].fPredictProbability = predict_np[i, 1]
        
    psm_new_list = re_rank(psm_rank_list)
    # del psm_list[:]
    # psm_list.extend(psm_new_list)

    psm_filtered_list_local = show_Fdr(psm_new_list, None, fdr=fdr_given)
    # show_Fdr_charge(psm_list)
    psm_filtered_list.extend(psm_filtered_list_local)
    
    return psm_filtered_list
    
    for i in range(4):
        sys.stdout.write('%.3f' % logreg.coef_[0][i])
        sys.stdout.write('\t')
    sys.stdout.write('\n')
        
    return psm_filtered_list

def test_voting(psm_list, fdr_given=None):
    # re_rank
    psm_new_list = []
    psm_dict = {}
    for oPsm in psm_list:
        oPsm.set_fdr_product()
        sId = oPsm.FileName + '_' + str(oPsm.ScanNumber)
        if sId in psm_dict:
            if oPsm.fdr_product < psm_dict[sId].fdr_product:
                psm_dict[sId] = oPsm
        else:
            psm_dict[sId] = oPsm
    
    for _key, value in psm_dict.items():
        psm_new_list.append(value)
    
    num_threshold = 20
    vote_threshold = 4
    fdr_threshold_list = [0.01] * num_threshold
    for i in range(num_threshold):
        fdr_threshold_list[i] += i * 0.001
    
    num_fwd_list = [0] * num_threshold
    num_rev_list = [0] * num_threshold
    num_shu_list = [0] * num_threshold
    for oPsm in psm_new_list:
        for i in range(num_threshold):
            if sum(val <= fdr_threshold_list[i] for val in oPsm.FDR_feature) >= vote_threshold:
                if oPsm.RealLabel == LabelFwd:
                    num_fwd_list[i] += 1
                elif oPsm.RealLabel == LabelRevTrain:
                    num_rev_list[i] += 1
                else:
                    num_shu_list[i] += 1
    print('\nFDR\t# FWD\t# REV\t#SHU')
    for i in range(num_threshold):
        print(str(fdr_threshold_list[i]) + '\t' + str(num_fwd_list[i]) + '\t' + str(num_rev_list[i]) + '\t' + str(num_shu_list[i]))
    
    return None

def writeout_feature(list_before_filtering, list_after_filtering):
    filename_before_str = '/media/xgo/Seagate/Proteomics/Experiments/BenchmarkRev/Sipros10/0820/fig2/fig2_Angelo_before.txt'
    filename_after_str = '/media/xgo/Seagate/Proteomics/Experiments/BenchmarkRev/Sipros10/0820/fig2/fig2_Angelo_after.txt'
    ms2_filename_dict = {}
    ms2_filename_str = ''
    index_int = 0
    scan_set = set()
    id_str = ''
    with open(filename_before_str, 'w') as fw:
        fw.write('RealLabel\tMVH\tXcorr\tWDP\tDiffMVH\tDiffXcorr\tDiffWDP\tMassError\tMissedCleavage\tOPSC\tSPSC\n')
        for oPsm in list_before_filtering:
            if oPsm.FileName in ms2_filename_dict:
                ms2_filename_str = ms2_filename_dict[oPsm.FileName]
            else:
                ms2_filename_dict[oPsm.FileName] = str(index_int)
                ms2_filename_str = str(index_int)
                index_int += 1
            id_str = ms2_filename_str + '_' + str(oPsm.ScanNumber) + '_' + str(oPsm.ParentCharge) + '_' + oPsm.IdentifiedPeptide
            if id_str in scan_set:
                continue
            else:
                scan_set.add(id_str)
                fw.write(str(oPsm.RealLabel))
                fw.write('\t')
                fw.write(str(oPsm.lfScores[0]))
                fw.write('\t')
                fw.write(str(oPsm.lfScores[1]))
                fw.write('\t')
                fw.write(str(oPsm.lfScores[2]))
                fw.write('\t')
                fw.write(str(oPsm.score_differential_list[9]))
                fw.write('\t')
                fw.write(str(oPsm.score_differential_list[10]))
                fw.write('\t')
                fw.write(str(oPsm.score_differential_list[11]))
                # fw.write('\t')
                # fw.write(str(abs(oPsm.iMassWindow)))
                fw.write('\t')
                fw.write(str(abs(oPsm.fMassDiff)))
                fw.write('\t')
                fw.write(str(oPsm.NMC))
                # fw.write('\t')
                # fw.write(str(oPsm.IPSC))
                fw.write('\t')
                fw.write(str(oPsm.OPSC))
                # fw.write('\t')
                # fw.write(str(oPsm.UPSC))
                fw.write('\t')
                fw.write(str(oPsm.SPSC))
                fw.write('\n')
    
    with open(filename_after_str, 'w') as fw:
        fw.write('RealLabel\tMVH\tXcorr\tWDP\tDiffMVH\tDiffXcorr\tDiffWDP\tMassError\tMissedCleavage\tOPSC\tSPSC\n')
        for oPsm in list_after_filtering:
            if oPsm.FileName in ms2_filename_dict:
                ms2_filename_str = ms2_filename_dict[oPsm.FileName]
            else:
                ms2_filename_dict[oPsm.FileName] = str(index_int)
                ms2_filename_str = str(index_int)
                index_int += 1
            fw.write(str(oPsm.RealLabel))
            fw.write('\t')
            fw.write(str(oPsm.lfScores[0]))
            fw.write('\t')
            fw.write(str(oPsm.lfScores[1]))
            fw.write('\t')
            fw.write(str(oPsm.lfScores[2]))
            fw.write('\t')
            fw.write(str(oPsm.score_differential_list[9]))
            fw.write('\t')
            fw.write(str(oPsm.score_differential_list[10]))
            fw.write('\t')
            fw.write(str(oPsm.score_differential_list[11]))
            fw.write('\t')
            fw.write(str(abs(oPsm.iMassWindow)))
            fw.write('\t')
            fw.write(str(abs(oPsm.fMassDiff)))
            fw.write('\t')
            fw.write(str(oPsm.NMC))
            # fw.write('\t')
            # fw.write(str(oPsm.IPSC))
            fw.write('\t')
            fw.write(str(oPsm.OPSC))
            # fw.write('\t')
            # fw.write(str(oPsm.UPSC))
            fw.write('\t')
            fw.write(str(oPsm.SPSC))
            fw.write('\n')

def logistic_regression_no_category(psm_dict, psm_list, psm_neg_list, fdr_given=None, psm_left_list=None, output_folder=None, input_file=None, config_dict=None):
    # machine learning
    # # construct training data
    #psm_list_selected = []
    train_data_list = []
    train_label_list = []
    test_data_list = []
    '''
    for key, lPsm in psm_dict.items():
        #sys.stdout.write(key + "\t")
        for oPsm in lPsm:
            if len(oPsm.feature_list) == 0:
                print 'check'
            if oPsm.TrainingLabel == LabelFiltered:
                continue
            test_data_list.append(oPsm.feature_list)
            if oPsm.TrainingLabel != LabelUnknown:
                train_data_list.append(oPsm.feature_list)
                train_label_list.append(oPsm.TrainingLabel)
    '''
    psm_rank_list = []
    '''
    psm_rank_U_list = []
    psm_rank_M_list = []
    psm_rank_D_list = []
    '''
    num_feature_int = (31 + len(ptm_selection_list))
    positive_int = 1
    negative_int = 0
    bDisableLocalRank = False
    if len(psm_list) < 800000:
        bDisableLocalRank = True
    for oPsm in psm_list:
        if len(oPsm.feature_list) != num_feature_int:
            pass
            # print 'check'
        if oPsm.TrainingLabel == LabelFiltered:
            continue
        if oPsm.RealLabel == LabelRevReserve:
            # train_data_list.append(oPsm.feature_list)
            # train_label_list.append(negative_int)
            continue
        # if oPsm.iLocalRank == 0 or oPsm.iLocalRank == 1:
        test_data_list.append(oPsm.feature_list)
        psm_rank_list.append(oPsm)
        '''
        if oPsm.iLocalRank == 0:
            psm_rank_U_list.append(oPsm)
        elif oPsm.iLocalRank == 2:
            psm_rank_D_list.append(oPsm)
        else:
            psm_rank_M_list.append(oPsm)
        '''
        # if oPsm.TrainingLabel != LabelUnknown and oPsm.iLocalRank == 0:
        # if oPsm.RealLabel == LabelFwd and (oPsm.iLocalRank == 0 ):
        if oPsm.RealLabel != LabelRevTrain and (oPsm.iLocalRank == 0 or bDisableLocalRank ):
            train_data_list.append(oPsm.feature_list)
            train_label_list.append(positive_int)
        elif oPsm.RealLabel == LabelRevTrain: # and (oPsm.iLocalRank == 0 or oPsm.iLocalRank == 1):
            train_data_list.append(oPsm.feature_list)
            train_label_list.append(negative_int)
                
    sys.stdout.write(str(len(train_data_list)) + "\t")
        
    train_data_np = np.array(train_data_list)[:, feature_selection_list]
    train_label_np = np.array(train_label_list)
    
    # only forward left
    unique_np = np.unique(train_label_np)
    if unique_np.shape[0] == 1:
        psm_filtered_list_local = show_Fdr(psm_list, None, fdr=fdr_given)
        return psm_filtered_list_local
    
    # # training
    # num_positive = float((train_label_np==LabelPositive).sum())
    # num_negative = float((train_label_np==LabelNegative).sum())
    # num_positive = 100.0
    # num_negative = 1.0
    # class_weight_dict = {0: (num_positive/(num_negative+num_positive)), 1:(num_negative/(num_negative+num_positive))}
    logreg = linear_model.LogisticRegression(penalty='l2', 
                                             dual=False, 
                                             C=1, 
                                             fit_intercept=True, 
                                             intercept_scaling=1, 
                                             # class_weight='balanced',
                                             # class_weight=class_weight_dict, 
                                             random_state=None, 
                                             solver='liblinear', 
                                             max_iter=100, 
                                             multi_class='ovr', 
                                             verbose=0, 
                                             warm_start=False, 
                                             n_jobs=-1)
    logreg.fit(train_data_np, train_label_np)

    for i in range(len(feature_selection_list)):
        sys.stdout.write('%.3f\n' % logreg.coef_[0][i])
    
    sys.stdout.write('\n')


    # # test
    test_unknown_np = np.array(test_data_list)[:, feature_selection_list]
    predict_np = logreg.predict_proba(test_unknown_np)
    
    for i in range(len(predict_np)):
        psm_rank_list[i].fPredictProbability = predict_np[i, 1]
        psm_rank_list[i].ML_feature.append(predict_np[i, 1])
    # fdr_rank(psm_list, 0)
    # return None
    '''
    # U
    test_data_list = []
    for oPsm in psm_rank_U_list:
        test_data_list.append(oPsm.feature_list)
    test_unknown_np = np.array(test_data_list)[:, feature_selection_list]
    predict_np = logreg.predict_proba(test_unknown_np)
    for i in range(len(predict_np)):
        psm_rank_U_list[i].fPredictProbability = predict_np[i, 1]
    psm_rank_U_list = re_rank(psm_rank_U_list)
    psm_rank_U_list = show_Fdr(psm_rank_U_list, None, fdr=fdr_given)
    # M
    test_data_list = []
    for oPsm in psm_rank_M_list:
        test_data_list.append(oPsm.feature_list)
    test_unknown_np = np.array(test_data_list)[:, feature_selection_list]
    predict_np = logreg.predict_proba(test_unknown_np)
    for i in range(len(predict_np)):
        psm_rank_M_list[i].fPredictProbability = predict_np[i, 1]
    psm_rank_M_list = re_rank(psm_rank_M_list)
    psm_rank_M_list = show_Fdr(psm_rank_M_list, None, fdr=fdr_given)
    # D
    test_data_list = []
    for oPsm in psm_rank_D_list:
        test_data_list.append(oPsm.feature_list)
    test_unknown_np = np.array(test_data_list)[:, feature_selection_list]
    predict_np = logreg.predict_proba(test_unknown_np)
    for i in range(len(predict_np)):
        psm_rank_D_list[i].fPredictProbability = predict_np[i, 1]
    psm_rank_D_list = re_rank(psm_rank_D_list)
    psm_rank_D_list = show_Fdr(psm_rank_D_list, None, fdr=fdr_given)
    psm_rank_list = psm_rank_U_list
    psm_rank_list.extend(psm_rank_M_list)
    psm_rank_list.extend(psm_rank_D_list)
    '''
    psm_new_list = re_rank(psm_rank_list)
    # del psm_list[:]
    # psm_list.extend(psm_new_list)
    
    # psm_new_list = mass_filter(psm_new_list)
    
    for fdr_f in [0.005, 0.01, 0.02]:
        temp_list = show_Fdr_varied(psm_new_list, fdr_f/3.0)
        '''
        folder_str = output_folder + 'fdr_' + str(fdr_f) +'/'
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
        generate_psm_pep_txt(input_file, folder_str, temp_list)
        '''
    for idx, num in enumerate(psm_num_list):
        sys.stdout.write('{:,d} '.format(num))
        sys.stdout.write('({:,d})\t'.format(true_psm_num_list[idx]))
    sys.stdout.write('\n')
    for idx, num in enumerate(pep_num_list):
        sys.stdout.write('{:,d} '.format(num))
        sys.stdout.write('({:,d})\t'.format(true_pep_num_list[idx]))
    sys.stdout.write('\n')
    
    
    if config_dict != None:
        if config_dict[FDR_Filtering_str] == 'PSM':
            psm_filtered_list_local = show_Fdr(psm_new_list, None, fdr=fdr_given)
        else:
            psm_filtered_list_local = show_Fdr_Pep(psm_new_list, None, fdr=fdr_given)
        # writeout_feature(psm_list, psm_filtered_list_local)
        return psm_filtered_list_local
    
    psm_filtered_list_local = show_Fdr(psm_new_list, None, fdr=fdr_given)
        # psm_filtered_list_local = show_Fdr_peptide(psm_list, None, fdr=fdr_given)
    # show_Fdr_charge(psm_list)
    
    
    
    psm_filtered_list = []
    psm_filtered_list.extend(psm_filtered_list_local)
    b = logreg.coef_[0]
    sum_float = 0
    for val in b:
        sum_float += abs(val)
    idx = 0
    for i in range(len(feature_name_list)):
        if i in feature_selection_list:
            # sys.stdout.write('%.3f' % logreg.coef_[0][idx])
            sys.stdout.write('%.5f' % (b[idx]/sum_float))
            idx += 1
        sys.stdout.write('\t')
    sys.stdout.write('\n')
    '''
    order_list = [2, 0, 1, 6, 4, 5, 12, 3, 7, 8, 9, 10, 11]
    print "\nfeature importance:"
    for ind in order_list:
        sys.stdout.write('%.5f\n' % logreg.coef_[0][ind])
    '''    
        
    '''
    iC = 0
    for oPsm in psm_filtered_list:
        if oPsm.ScoreAgreement == 1 and oPsm.RealLabel == LabelFwd:
            iC += 1
    print str(iC)
    '''        
        
    return psm_filtered_list

def logistic_regression_1_LR(psm_dict, psm_list, psm_neg_list, fdr_given=None, psm_left_list=None, output_folder=None, input_file=None):
    # machine learning
    # # construct training data
    #psm_list_selected = []
    data_list = []
    label_list = []
    unknown_list = []
    id_list = []
    
    positive_int = 1
    negative_int = 0
    
    test_list_1 = []
    train_list_1 = []
    label_list_1 = []
    psm_list_1 = []
    
    test_list_2 = []
    train_list_2 = []
    label_list_2 = []
    psm_list_2 = []
    
    psm_list_3 = []
    test_list_3 = []
    
    num_feature_int = (31 + len(ptm_selection_list))           
    for oPsm in psm_list:
        if len(oPsm.feature_list) != num_feature_int:
            print('check')
        if oPsm.iLocalRank == 0: # U
            test_list_1.append(oPsm.feature_list)
            psm_list_1.append(oPsm)
            train_list_1.append(oPsm.feature_list)
            if oPsm.RealLabel != LabelRevTrain:
                label_list_1.append(positive_int)
                # train_list_2.append(oPsm.feature_list)
                # label_list_2.append(positive_int)
            else:
                label_list_1.append(negative_int)
        else:
            if oPsm.iLocalRank == 1:
                test_list_2.append(oPsm.feature_list)
                psm_list_2.append(oPsm)
                if oPsm.RealLabel != LabelRevTrain:
                    # pass
                    train_list_2.append(oPsm.feature_list)
                    label_list_2.append(positive_int)
                else:
                    train_list_2.append(oPsm.feature_list)
                    label_list_2.append(negative_int)
            elif oPsm.iLocalRank == 3:
                test_list_2.append(oPsm.feature_list)
                psm_list_2.append(oPsm)
                if oPsm.RealLabel != LabelRevTrain:
                    pass
                else:
                    train_list_2.append(oPsm.feature_list)
                    label_list_2.append(negative_int)
            else:
                test_list_3.append(oPsm.feature_list)
                psm_list_3.append(oPsm)
                if oPsm.RealLabel != LabelRevTrain:
                    pass
                else:
                    pass
                    # train_list_1.append(oPsm.feature_list)
                    # label_list_1.append(negative_int)
    
    # # training
    print('train 1 size: %d' % len(train_list_1))
    train_data_np = np.array(train_list_1)[:, feature_selection_list]
    train_label_np = np.array(label_list_1)
    logreg = linear_model.LogisticRegression(penalty='l2', 
                                             dual=False, 
                                             C=1, 
                                             fit_intercept=True, 
                                             intercept_scaling=1, 
                                             # class_weight='balanced',
                                             # class_weight=class_weight_dict, 
                                             random_state=None, 
                                             solver='liblinear', 
                                             max_iter=100, 
                                             multi_class='ovr', 
                                             verbose=0, 
                                             warm_start=False, 
                                             n_jobs=-1)
    logreg.fit(train_data_np, train_label_np)
    
    # # test
    test_data_np = np.array(test_list_1)[:, feature_selection_list]
    predict_np = logreg.predict_proba(test_data_np)
    for i in range(len(predict_np)):
        psm_list_1[i].fPredictProbability = predict_np[i, 1]
        
    psm_new_list = re_rank(psm_list_1)
    psm_list_1 = show_Fdr(psm_new_list, None, fdr=fdr_given)
    
    test_data_np = np.array(test_list_3)[:, feature_selection_list]
    predict_np = logreg.predict_proba(test_data_np)
    for i in range(len(predict_np)):
        psm_list_3[i].fPredictProbability = predict_np[i, 1]
        
    psm_new_list = re_rank(psm_list_3)
    psm_list_3 = show_Fdr(psm_new_list, None, fdr=fdr_given)  
        
    # # training
    print('train 2 size: %d' % len(train_list_2))
    train_data_np = np.array(train_list_2)[:, feature_selection_list]
    train_label_np = np.array(label_list_2)
    logreg = linear_model.LogisticRegression(penalty='l2', 
                                             dual=False, 
                                             C=1, 
                                             fit_intercept=True, 
                                             intercept_scaling=1, 
                                             # class_weight='balanced',
                                             # class_weight=class_weight_dict, 
                                             random_state=None, 
                                             solver='liblinear', 
                                             max_iter=100, 
                                             multi_class='ovr', 
                                             verbose=0, 
                                             warm_start=False, 
                                             n_jobs=-1)
    logreg.fit(train_data_np, train_label_np)
    
    # # test
    test_data_np = np.array(test_list_2)[:, feature_selection_list]
    predict_np = logreg.predict_proba(test_data_np)
    for i in range(len(predict_np)):
        psm_list_2[i].fPredictProbability = predict_np[i, 1]
        
    psm_new_list = re_rank(psm_list_2)
    psm_list_2 = show_Fdr(psm_new_list, None, fdr=fdr_given)  
    
    psm_list_all = []
    psm_list_all.extend(psm_list_1)
    psm_list_all.extend(psm_list_2)
    psm_list_all.extend(psm_list_3)
    psm_new_list = re_rank(psm_list_all)
    psm_list_all = show_Fdr(psm_list_all, None, fdr=fdr_given)


def logistic_regression_2_LR(psm_dict, psm_list, psm_neg_list, fdr_given=None, psm_left_list=None, output_folder=None, input_file=None):
    # machine learning
    # # construct training data
    #psm_list_selected = []
    data_list = []
    label_list = []
    unknown_list = []
    id_list = []
    
    positive_int = 1
    negative_int = 0
    
    test_list_1 = []
    train_list_1 = []
    label_list_1 = []
    psm_list_1 = []
    
    test_list_2 = []
    train_list_2 = []
    label_list_2 = []
    psm_list_2 = []
    
    psm_list_3 = []
    test_list_3 = []
    
    num_feature_int = (31 + len(ptm_selection_list))           
    for oPsm in psm_list:
        if len(oPsm.feature_list) != num_feature_int:
            # pass
            print('check')
        if oPsm.iLocalRank == 0: # U
            test_list_1.append(oPsm.feature_list)
            psm_list_1.append(oPsm)
            train_list_1.append(oPsm.feature_list)
            if oPsm.RealLabel != LabelRevTrain:
                label_list_1.append(positive_int)
                # train_list_2.append(oPsm.feature_list)
                # label_list_2.append(positive_int)
            else:
                label_list_1.append(negative_int)
        else:
            if oPsm.iLocalRank == 1:
                test_list_2.append(oPsm.feature_list)
                psm_list_2.append(oPsm)
                if oPsm.RealLabel != LabelRevTrain:
                    # pass
                    train_list_2.append(oPsm.feature_list)
                    label_list_2.append(positive_int)
                else:
                    train_list_2.append(oPsm.feature_list)
                    label_list_2.append(negative_int)
            elif oPsm.iLocalRank == 3:
                test_list_2.append(oPsm.feature_list)
                psm_list_2.append(oPsm)
                if oPsm.RealLabel != LabelRevTrain:
                    pass
                else:
                    train_list_2.append(oPsm.feature_list)
                    label_list_2.append(negative_int)
            else:
                test_list_3.append(oPsm.feature_list)
                psm_list_3.append(oPsm)
                if oPsm.RealLabel != LabelRevTrain:
                    pass
                else:
                    pass
                    # train_list_1.append(oPsm.feature_list)
                    # label_list_1.append(negative_int)
    
    # # training
    print('train 1 size: %d' % len(train_list_1))
    train_data_np = np.array(train_list_1)[:, feature_selection_list]
    train_label_np = np.array(label_list_1)
    logreg = linear_model.LogisticRegression(penalty='l2', 
                                             dual=False, 
                                             C=1, 
                                             fit_intercept=True, 
                                             intercept_scaling=1, 
                                             # class_weight='balanced',
                                             # class_weight=class_weight_dict, 
                                             random_state=None, 
                                             solver='liblinear', 
                                             max_iter=100, 
                                             multi_class='ovr', 
                                             verbose=0, 
                                             warm_start=False, 
                                             n_jobs=-1)
    logreg.fit(train_data_np, train_label_np)
    
    # # test
    test_data_np = np.array(test_list_1)[:, feature_selection_list]
    predict_np = logreg.predict_proba(test_data_np)
    for i in range(len(predict_np)):
        psm_list_1[i].fPredictProbability = predict_np[i, 1]
        
    psm_new_list = re_rank(psm_list_1)
    psm_list_1 = show_Fdr(psm_new_list, None, fdr=fdr_given)
    
    test_data_np = np.array(test_list_3)[:, feature_selection_list]
    predict_np = logreg.predict_proba(test_data_np)
    for i in range(len(predict_np)):
        psm_list_3[i].fPredictProbability = predict_np[i, 1]
        
    psm_new_list = re_rank(psm_list_3)
    psm_list_3 = show_Fdr(psm_new_list, None, fdr=fdr_given)  
        
    # # training
    print('train 2 size: %d' % len(train_list_2))
    train_data_np = np.array(train_list_2)[:, feature_selection_list]
    train_label_np = np.array(label_list_2)
    logreg = linear_model.LogisticRegression(penalty='l2', 
                                             dual=False, 
                                             C=1, 
                                             fit_intercept=True, 
                                             intercept_scaling=1, 
                                             # class_weight='balanced',
                                             # class_weight=class_weight_dict, 
                                             random_state=None, 
                                             solver='liblinear', 
                                             max_iter=100, 
                                             multi_class='ovr', 
                                             verbose=0, 
                                             warm_start=False, 
                                             n_jobs=-1)
    logreg.fit(train_data_np, train_label_np)
    
    # # test
    test_data_np = np.array(test_list_2)[:, feature_selection_list]
    predict_np = logreg.predict_proba(test_data_np)
    for i in range(len(predict_np)):
        psm_list_2[i].fPredictProbability = predict_np[i, 1]
        
    psm_new_list = re_rank(psm_list_2)
    psm_list_2 = show_Fdr(psm_new_list, None, fdr=fdr_given)  
    
    psm_list_all = []
    psm_list_all.extend(psm_list_1)
    psm_list_all.extend(psm_list_2)
    psm_list_all.extend(psm_list_3)
    psm_new_list = re_rank(psm_list_all)
    psm_list_all = show_Fdr(psm_list_all, None, fdr=fdr_given)

import random

def test_train_bias(psm_dict, psm_list, psm_neg_list, fdr_given=None, psm_left_list=None):
    # machine learning
    # # construct training data
    train_data_list = []
    train_label_list = []
    train_id_list = []
    test_data_list = []
    test_id_list = []
    # # psm with the same name
    psm_sa_1_dict = {}
    
    for oPsm in psm_list:
        if oPsm.ScoreAgreement <= 1:
            sId = oPsm.FileName + '_' + str(oPsm.ScanNumber)
            if sId in psm_sa_1_dict:
                if psm_sa_1_dict[sId] == True:
                    train_data_list.append(oPsm.feature_list)
                    if oPsm.RealLabel == LabelRevTrain:
                        train_label_list.append(0)
                    else:
                        train_label_list.append(1)
                    train_id_list.append(oPsm.iInnerId)
                else:
                    test_data_list.append(oPsm.feature_list)
                    test_id_list.append(oPsm.iInnerId)
            else:
                f_rand = random.random()
                if f_rand >= 0.5:
                    psm_sa_1_dict[sId] = True
                    train_data_list.append(oPsm.feature_list)
                    if oPsm.RealLabel == LabelRevTrain:
                        train_label_list.append(0)
                    else:
                        train_label_list.append(1)
                    train_id_list.append(oPsm.iInnerId)
                else:
                    psm_sa_1_dict[sId] = False
                    test_data_list.append(oPsm.feature_list)
                    test_id_list.append(oPsm.iInnerId)
        else:
            f_rand = random.random()
            if f_rand >= 0.5:
                train_data_list.append(oPsm.feature_list)
                if oPsm.RealLabel == LabelRevTrain:
                    train_label_list.append(0)
                else:
                    train_label_list.append(1)
                train_id_list.append(oPsm.iInnerId)
            else:
                test_data_list.append(oPsm.feature_list)
                test_id_list.append(oPsm.iInnerId)
                
    train_data_np = np.array(train_data_list)
    sys.stdout.write(str(len(train_data_list)) + "\t")
    train_label_np = np.array(train_label_list)
        
    train_data_np = train_data_np[:, feature_selection_list]
   
    logreg = linear_model.LogisticRegression(penalty='l2', 
                                             dual=False, 
                                             C=1, 
                                             fit_intercept=True, 
                                             intercept_scaling=1, 
                                             # class_weight='balanced',
                                             # class_weight=class_weight_dict, 
                                             random_state=None, 
                                             solver='liblinear', 
                                             max_iter=100, 
                                             multi_class='ovr', 
                                             verbose=0, 
                                             warm_start=False, 
                                             n_jobs=-1)
    logreg.fit(train_data_np, train_label_np)

    # # test
    test_data_np = np.array(test_data_list)
    test_unknown_np = test_data_np[:, feature_selection_list]
    predict_np = logreg.predict_proba(test_unknown_np)
    #predict_np = logreg.predict_proba(train_data_np)
    for i in range(len(predict_np)):
        psm_list[test_id_list[i]].fPredictProbability = predict_np[i, 1]
        
    psm_new_list = re_rank(psm_list)
    print('testing results:')
    for i in range(1, 11):
        fdr_f = 0.001 * i
        show_Fdr_varied(psm_new_list, fdr_f)
    
    for i in range(len(predict_np)):
        psm_list[test_id_list[i]].fPredictProbability = 0 
    predict_np = logreg.predict_proba(train_data_np)
    for i in range(len(predict_np)):
        psm_list[train_id_list[i]].fPredictProbability = predict_np[i, 1]
    print('train results:')
    psm_new_list = re_rank(psm_list)
    for i in range(1, 11):
        fdr_f = 0.001 * i
        show_Fdr_varied(psm_new_list, fdr_f)       
        
    return None

from sklearn import svm

def svm_one_class_test(psm_dict, psm_list, psm_neg_list, fdr_given=None, psm_left_list=None):
    # machine learning
    # # construct training data
    psm_filtered_list = []
    #psm_list_selected = []
    data_list = []
    label_list = []
    unknown_list = []
    id_list = [] 
    for oPsm in psm_list:
        if len(oPsm.feature_list) == 0:
            print('check')
        if oPsm.TrainingLabel == LabelFiltered:
            continue
        unknown_list.append(oPsm.feature_list)
        id_list.append(oPsm.iInnerId)
        if oPsm.TrainingLabel == LabelNegative:
            data_list.append(oPsm.feature_list)
            label_list.append(oPsm.TrainingLabel)
    
           
    for oPsm in psm_neg_list:
        data_list.append(oPsm.feature_list)
        label_list.append(oPsm.TrainingLabel)
                
    data_np = np.array(data_list)
    sys.stdout.write(str(len(data_list)) + "\t")
    label_np = np.array(label_list)
    train_data_np = data_np[:, feature_selection_list]
    
    clf = svm.OneClassSVM(nu=0.1, kernel="rbf", gamma=0.1)
    
    clf.fit(train_data_np)
    #predict_np = logreg.predict(train_data_np)

    # # test
    unknown_np = np.array(unknown_list)
    test_unknown_np = unknown_np[:, feature_selection_list]
    predict_np = clf.predict_proba(test_unknown_np)
    #predict_np = logreg.predict_proba(train_data_np)
    for i in range(len(predict_np)):
        psm_list[id_list[i]].fPredictProbability = predict_np[i, 1]
        
    psm_new_list = re_rank(psm_list)
    del psm_list[:]
    psm_list.extend(psm_new_list)
    # np.savetxt("predict_probability.txt", predict_np)
        #print str(len(psm_list_selected))
    if not psm_left_list is None:
        psm_filtered_list_local = show_Fdr(psm_left_list, None, fdr=fdr_given)
    else:
        psm_filtered_list_local = show_Fdr(psm_list, None, fdr=fdr_given)
    # show_Fdr_charge(psm_list)
    psm_filtered_list.extend(psm_filtered_list_local)
    sys.stdout.write('\n')
        
    return psm_filtered_list

def logistic_regression_no_category_rt(psm_dict, psm_list, psm_neg_list, fdr_given=None, psm_left_list=None):
    # machine learning
    # # construct training data
    psm_filtered_list = []
    #psm_list_selected = []
    data_list = []
    label_list = []
    unknown_list = []
    id_list = []
    for key, lPsm in psm_dict.items():
        #sys.stdout.write(key + "\t")
        for oPsm in lPsm:
            if len(oPsm.feature_list) == 0:
                print('check')
                
            if oPsm.RealLabel == LabelRevTrain:
                data_list.append(oPsm.feature_list)
                label_list.append(oPsm.TrainingLabel)
            if oPsm.TrainingLabel == LabelFiltered:
                if oPsm.RealLabel != LabelRevTrain:
                    data_list.append(oPsm.feature_list)
                    label_list.append(oPsm.TrainingLabel)
                continue
            unknown_list.append(oPsm.feature_list)
            id_list.append(oPsm.iInnerId)
           
    for oPsm in psm_neg_list:
        data_list.append(oPsm.feature_list)
        label_list.append(oPsm.TrainingLabel)
                
    data_np = np.array(data_list)
    sys.stdout.write(str(len(data_list)) + "\t")
    label_np = np.array(label_list)
        
    train_data_np = data_np[:, feature_selection_list]
    '''
    if not psm_left_list is None:
        np.savetxt("/media/xgo/Seagate/Proteomics/Experiments/Angelo/Sipros10/RTime/check/train_data_np.txt", train_data_np)
    '''
    # # training
    #class_dict = {0: 1, 1:1000}
    logreg = linear_model.LogisticRegression(penalty='l2', 
                                             dual=False, 
                                             C=1, 
                                             fit_intercept=True, 
                                             intercept_scaling=1, 
                                             class_weight='balanced', 
                                             random_state=None, 
                                             solver='liblinear', 
                                             max_iter=100, 
                                             multi_class='ovr', 
                                             verbose=0, 
                                             warm_start=False, 
                                             n_jobs=-1)
    logreg.fit(train_data_np, label_np)
    #predict_np = logreg.predict(train_data_np)

    # # test
    unknown_np = np.array(unknown_list)
    test_unknown_np = unknown_np[:, feature_selection_list]
    predict_np = logreg.predict_proba(test_unknown_np)
    #predict_np = logreg.predict_proba(train_data_np)
    for i in range(len(predict_np)):
        psm_list[id_list[i]].fPredictProbability = predict_np[i, 1]

    # np.savetxt("predict_probability.txt", predict_np)
        #print str(len(psm_list_selected))
    if not psm_left_list is None:
        psm_filtered_list_local = show_Fdr(psm_left_list, key, fdr=fdr_given)
    else:
        psm_filtered_list_local = show_Fdr(psm_list, key, fdr=fdr_given)
    psm_filtered_list.extend(psm_filtered_list_local)
    # sys.stdout.write('\n')
    # show_Fdr_category(psm_dict)
    #print 'Coefficient of the features in the decision function:'
    #print logreg.coef_
    idx = 0
    for i in range(len(feature_name_list)):
        if i in feature_selection_list:
            sys.stdout.write('%.3f' % logreg.coef_[0][idx])
            idx += 1
        sys.stdout.write('\t')
    sys.stdout.write('\n')
        
    return psm_filtered_list

def get_num_missed_cleavage_sites(sIdentifiedSeq, sResiduesBeforeCleavage, sResiduesAfterCleavage):
    count_int = 0
    for i in range(len(sIdentifiedSeq) - 1):
        if sIdentifiedSeq[i] in sResiduesBeforeCleavage and sIdentifiedSeq[i + 1] in sResiduesAfterCleavage:
            count_int += 1
    return count_int

def generate_pep_pro_features(lPsm, config_dict):
    # number of replicate spectra
    identified_peptide_dict = {}
    # number of sibling ions
    identified_peptide_with_different_charges_dict = {}
    # peptide to protein dictionary
    original_pep_pro_dict = {}
    for oPsm in lPsm:
        # number of missed cleavage sites
        '''
        oPsm.NMC = get_num_missed_cleavage_sites(oPsm.OriginalPeptide, 
                                                 config_dict[cleave_after_residues_str], 
                                                 config_dict[cleave_before_residues_str])
        '''
        if oPsm.IdentifiedPeptide in identified_peptide_with_different_charges_dict:
            if oPsm.ParentCharge not in identified_peptide_with_different_charges_dict[oPsm.IdentifiedPeptide]:
                identified_peptide_with_different_charges_dict[oPsm.IdentifiedPeptide].append(oPsm.ParentCharge)
        else:
            identified_peptide_with_different_charges_dict[oPsm.IdentifiedPeptide] = [oPsm.ParentCharge]
        if oPsm.OriginalPeptide not in original_pep_pro_dict:
            original_pep_pro_dict[oPsm.OriginalPeptide] = set()
            for pro_str in oPsm.protein_list:
                original_pep_pro_dict[oPsm.OriginalPeptide].add(pro_str)
        else:
            for pro_str in oPsm.protein_list:
                original_pep_pro_dict[oPsm.OriginalPeptide].add(pro_str)
        sId = oPsm.IdentifiedPeptide + '_' + str(oPsm.ParentCharge)
        if sId in identified_peptide_dict:
            identified_peptide_dict[sId] += 1
        else:
            identified_peptide_dict[sId] = 1
        
    # number of sibling modification
    original_peptide_dict = {}
    pattern = re.compile('[\W_]+')
    for (pep_str, charge_list) in identified_peptide_with_different_charges_dict.items():
        # get the original peptide
        pep_original_str = '[' + pattern.sub('', pep_str) + ']'
        if pep_original_str not in original_peptide_dict:
            original_peptide_dict[pep_original_str] = 1
        else:
            original_peptide_dict[pep_original_str] += 1
    # number of sibling ions
    for oPsm in lPsm:
        if oPsm.IdentifiedPeptide not in identified_peptide_with_different_charges_dict:
            print('error in number of sibling ions')
        else:
            pass
            # oPsm.IPSC = len(identified_peptide_with_different_charges_dict[oPsm.IdentifiedPeptide])
        # # get the number of sibling modification
        if oPsm.OriginalPeptide not in original_peptide_dict:
            print('error in number of sibling modification')
        else:
            pass
            # oPsm.OPSC = original_peptide_dict[oPsm.OriginalPeptide]
        # # get the number of replicate spectra
        sId = oPsm.IdentifiedPeptide + '_' + str(oPsm.ParentCharge)
        if sId in identified_peptide_dict:
            oPsm.NRS = identified_peptide_dict[sId]
        else:
            print('error in number of replicate spectra')
    
    # protein to peptide dictionary
    pro_original_pep_dict = {}
    for (pep_str, pro_list) in original_pep_pro_dict.items():
        for pro_str in pro_list:
            if pro_str not in pro_original_pep_dict:
                pro_original_pep_dict[pro_str] = set()
                pro_original_pep_dict[pro_str].add(pep_str)
            else:
                pro_original_pep_dict[pro_str].add(pep_str)
    
    # number of sibling unique peptide and number of sibling shared peptide
    '''
    pro_unique_dict = {}
    pro_shared_dict = {}
    changed_flag = False
    for oPsm in lPsm:
        pro_list = original_pep_pro_dict[oPsm.OriginalPeptide]
        changed_flag = False
        for protein in pro_list:
            if not protein in oPsm.protein_list:
                oPsm.protein_list.append(protein)
                changed_flag = True
        if changed_flag:
            oPsm.set_protein_names()
            oPsm.RealLabel = protein_type(oPsm.ProteinNames)
        if len(pro_list) > 1:
            for protein in pro_list:
                if protein in pro_shared_dict: 
                    pro_shared_dict[protein] += 1
                else:
                    pro_shared_dict[protein] = 1
        else:
            pro_str = pro_list.pop()
            pro_list.add(pro_str)
            if pro_str in pro_unique_dict:
                pro_unique_dict[pro_str] += 1
            else:
                pro_unique_dict[pro_str] = 1
    
    # collect features
    for oPsm in lPsm:
        pro_list = original_pep_pro_dict[oPsm.OriginalPeptide]
        for protein in pro_list:
            if protein in pro_unique_dict:
                oPsm.UPSC += pro_unique_dict[protein]
                if len(pro_list) == 1:
                    oPsm.UPSC -= 1
            if protein in pro_shared_dict:
                oPsm.SPSC += pro_shared_dict[protein]
                if len(pro_list) > 1:
                    oPsm.SPSC -= 1
    
    '''
    
    # prophet way
    num_unique_per_pro = 0
    num_shared_per_pro = 0
    num_balanced_shared_per_pro = 0
    
    max_unique_per_psm = 0
    max_shared_per_psm = 0
    max_balanced_shared_per_psm = 0
    
    num_unique_per_psm = 0
    num_shared_per_psm = 0
    num_balanced_shared_per_psm = 0
    
    num_per_pro = 0 # linked to a protein
    num_balanced_per_pro = 0 # linked to a protein
    
    max_per_psm = 0 # linked to protein
    max_balanced_per_psm = 0 # linked to a protein
    
    max_linked_unique_per_psm = 0
    max_linked_shared_per_psm = 0
    max_linked_balanced_unique_per_psm = 0
    max_linked_balanced_shared_per_psm = 0
    
    
    for oPsm in lPsm:
        num_unique_per_psm = 0
        num_shared_per_psm = 0
        num_balanced_shared_per_psm = 0
        max_per_psm = 0
        max_balanced_per_psm = 0
        max_unique_per_psm = 0
        max_shared_per_psm = 0
        max_balanced_shared_per_psm = 0
        # get the protein list first
        if oPsm.OriginalPeptide not in original_pep_pro_dict:
            print('error peptide not in the dictionary')
        else:
            pro_list = original_pep_pro_dict[oPsm.OriginalPeptide]
            for pro_str in pro_list:
                # debug
                num_unique_per_pro = 0
                num_shared_per_pro = 0
                num_balanced_shared_per_pro = 0
                num_per_pro = 0
                num_balanced_per_pro = 0
                # debug
                if pro_str not in pro_original_pep_dict:
                    print('error in number of sibling peptide')
                else:
                    for pep_str in pro_original_pep_dict[pro_str]:
                        if pep_str not in original_pep_pro_dict:
                            print('error in pep dictionary')
                        else:
                            if len(original_pep_pro_dict[pep_str]) == 1: # unique
                                num_unique_per_pro += 1
                            else: # shared
                                num_shared_per_pro += 1
                                num_balanced_shared_per_pro += 1.0/float(len(original_pep_pro_dict[pep_str]))
                    num_per_pro += num_unique_per_pro + num_shared_per_pro
                    num_balanced_per_pro += num_unique_per_pro + num_balanced_shared_per_pro
                num_unique_per_psm += num_unique_per_pro
                num_shared_per_psm += num_shared_per_pro
                num_balanced_shared_per_psm += num_balanced_shared_per_pro
                if num_per_pro > max_per_psm:
                    max_per_psm = num_per_pro
                    max_linked_unique_per_psm = num_unique_per_pro
                    max_linked_shared_per_psm = num_shared_per_pro
                if num_balanced_per_pro > max_balanced_per_psm:
                    max_balanced_per_psm = num_balanced_per_pro
                    max_linked_balanced_unique_per_psm = num_unique_per_pro
                    max_linked_balanced_shared_per_psm = num_balanced_shared_per_pro
                if num_unique_per_pro > max_unique_per_psm:
                    max_unique_per_psm = num_unique_per_pro
                if num_shared_per_pro > max_shared_per_psm:
                    max_shared_per_psm = num_shared_per_pro
                if num_balanced_shared_per_pro > max_balanced_shared_per_psm:
                    max_balanced_shared_per_psm = num_balanced_shared_per_pro
            
            if len(original_pep_pro_dict[oPsm.OriginalPeptide]) == 1:
                num_unique_per_psm -= 1
            else:
                num_shared_per_psm -= 1
        
        oPsm.UPSC = num_unique_per_psm
        oPsm.SPSC = num_shared_per_psm
        oPsm.SPSC = num_shared_per_psm + num_unique_per_psm
        oPsm.UPSC = max_unique_per_psm
        oPsm.SPSC = max_shared_per_psm
        oPsm.SPSC = max_shared_per_psm + max_unique_per_psm
        oPsm.SPSC = max_per_psm
        oPsm.UPSC = max_linked_unique_per_psm
        oPsm.SPSC = max_linked_shared_per_psm
        oPsm.SPSC = max_linked_shared_per_psm + max_linked_unique_per_psm
        # balanced
        # oPsm.UPSC = num_unique_per_psm
        # oPsm.SPSC = num_balanced_shared_per_psm
        # oPsm.SPSC = num_balanced_shared_per_psm + num_unique_per_psm
        oPsm.UPSC = max_linked_balanced_unique_per_psm
        oPsm.SPSC = max_linked_balanced_shared_per_psm
        oPsm.SPSC = max_linked_balanced_unique_per_psm + max_linked_balanced_shared_per_psm
        oPsm.SPSC = max_balanced_per_psm
        oPsm.SPSC = max_balanced_shared_per_psm + max_unique_per_psm
    
    

# generate OPSC (# sibling modification), IPSC (# sibling ions, charge), UPSC (# sibling peptides), NMC (# missed cleavage sites)
def generate_Prophet_features(lPsm, config_dict):
    # peptide with PTM dictionary is for IPSC
    peptide_with_modification_dict = {}
    # peptide without PTM dictionary is for OPSC
    peptide_dict = {}
    peptide_protein_dict = {}
    for oPsm in lPsm:
        oPsm.NMC = get_num_missed_cleavage_sites(oPsm.OriginalPeptide, 
                                                 config_dict[cleave_after_residues_str],
                                                 config_dict[cleave_before_residues_str])
        if oPsm.IdentifiedPeptide in peptide_with_modification_dict:
            peptide_with_modification_dict[oPsm.IdentifiedPeptide] += 1
        else:
            peptide_with_modification_dict[oPsm.IdentifiedPeptide] = 1
        if oPsm.OriginalPeptide in peptide_protein_dict:
            pro_list = peptide_protein_dict[oPsm.OriginalPeptide]
            for protein in oPsm.protein_list:
                if not protein in pro_list:
                    pro_list.append(protein)
        else:
            pro_list = []
            pro_list.extend(oPsm.protein_list)
            peptide_protein_dict[oPsm.OriginalPeptide] = pro_list
            

    pattern = re.compile('[^\w\[\]]')
    for key, _value in peptide_with_modification_dict.items():
        peptide_str = pattern.sub('', key)
        if peptide_str in peptide_dict:
            peptide_dict[peptide_str] += 1
        else:
            peptide_dict[peptide_str] = 1
    
    # # sibling peptides
    pro_unique_dict = {}
    pro_shared_dict = {}
    # debug
    pro_balanced_shared_dict = {}
    # debug
    changed_flag = False
    for oPsm in lPsm:
        pro_list = peptide_protein_dict[oPsm.OriginalPeptide]
        changed_flag = False
        for protein in pro_list:
            if not protein in oPsm.protein_list:
                oPsm.protein_list.append(protein)
                changed_flag = True
        if changed_flag:
            oPsm.set_protein_names()
            oPsm.RealLabel = protein_type(oPsm.ProteinNames)
            oPsm.ecoli_label = protein_type_ecoli(oPsm.protein_list)
        if len(pro_list) > 1:
            num_pro_f = float(len(pro_list))
            for protein in pro_list:
                if protein in pro_shared_dict: 
                    pro_shared_dict[protein] += 1
                else:
                    pro_shared_dict[protein] = 1
                # debug
                if protein in pro_balanced_shared_dict:
                    pro_balanced_shared_dict[protein] += 1.0/num_pro_f
                else:
                    pro_balanced_shared_dict[protein] = 1.0/num_pro_f
                # debug
        else:
            if pro_list[0] in pro_unique_dict:
                pro_unique_dict[pro_list[0]] += 1
            else:
                pro_unique_dict[pro_list[0]] = 1
    
    # collect features
    for oPsm in lPsm:
        oPsm.IPSC = peptide_with_modification_dict[oPsm.IdentifiedPeptide] - 1
        oPsm.OPSC = peptide_dict[oPsm.OriginalPeptide] - 1
        pro_list = peptide_protein_dict[oPsm.OriginalPeptide]
        for protein in pro_list:
            if protein in pro_unique_dict:
                oPsm.UPSC += pro_unique_dict[protein]
                if len(pro_list) == 1:
                    oPsm.UPSC -= 1
                
            if protein in pro_shared_dict:
                oPsm.SPSC += pro_shared_dict[protein]
                if len(pro_list) > 1:
                    oPsm.SPSC -= 1
                    
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
    # psm_set = set()
    for oPsm in lPsm:
        oPsm.NMC = get_num_missed_cleavage_sites(oPsm.OriginalPeptide, 
                                                 config_dict[cleave_after_residues_str],
                                                 config_dict[cleave_before_residues_str])
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
                if not protein in pro_list:
                    pro_list.append(protein)
        else:
            pro_list = []
            for pro in oPsm.protein_list:
                if pro not in pro_list:
                    pro_list.extend(oPsm.protein_list)
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
        changed_flag = False
        for protein in pro_list:
            if not protein in oPsm.protein_list:
                oPsm.protein_list.append(protein)
                changed_flag = True
        if len(oPsm.protein_list) != len(pro_list):
            print('check 3')
        
        if changed_flag:
            oPsm.set_protein_names()
            oPsm.RealLabel = protein_type(oPsm.ProteinNames)
            oPsm.ecoli_label = protein_type_ecoli(oPsm.protein_list)
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
        oPsm.SPSC = max_linked_unique_per_psm #+ max_linked_shared_per_psm
        # oPsm.SPSC = float(max_linked_unique_per_psm)/float(len(oPsm.protein_list)) + float(max_linked_shared_per_psm)/float(len(oPsm.protein_list))
        # oPsm.SPSC = max_per_psm
        
        # oPsm.UPSC = max_unique_per_psm
        # oPsm.SPSC = max_shared_per_psm
        
        '''
        # # balanced
        oPsm.UPSC = num_unique_per_psm
        oPsm.SPSC = num_balanced_shared_per_psm
        oPsm.SPSC = num_unique_per_psm + num_balanced_shared_per_psm
        oPsm.UPSC = max_unique_per_psm
        oPsm.SPSC = max_balanced_shared_per_psm
        oPsm.UPSC = max_linked_balanced_unique_per_psm
        oPsm.SPSC = max_linked_balanced_shared_per_psm
        # oPsm.SPSC = max_unique_per_psm + max_balanced_shared_per_psm
        # oPsm.SPSC = max_linked_balanced_unique_per_psm + max_linked_balanced_shared_per_psm
        '''
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
min_peptide_per_protein = 2
min_unique_peptide_per_protein = 1
remove_decoy_identification = 'No'

pep_iden_str = '[Peptide_Identification]'
fasta_database_str = 'FASTA_Database'
pro_iden_str = '[Protein_Identification]'
decoy_prefix_str = 'Decoy_Prefix'
FDR_Filtering_str = 'FDR_Filtering'
min_peptide_per_protein_str = 'Min_Peptide_Per_Protein'
min_unique_peptide_per_protein_str = 'Min_Unique_Peptide_Per_Protein'
remove_decoy_identification_str = 'Remove_Decoy_Identification'
cleave_after_residues_str = 'Cleave_After_Residues'
cleave_before_residues_str = 'Cleave_Before_Residues'

## Parse config file
def parse_config(config_filename):

    # Save config values to dictionary
    config_dict = {}    # initialize dictionay

    # Call Yinfeng's parseconfig.py module
    check_file_exist(config_filename)
    # Save all config values to dictionary
    all_config_dict = parseconfig.parseConfigKeyValues(config_filename)

    # only save protein_identification config info to config_dict
    config_dict[decoy_prefix_str] = decoy_prefix
    config_dict[min_peptide_per_protein_str] = min_peptide_per_protein
    config_dict[min_unique_peptide_per_protein_str] = min_unique_peptide_per_protein
    config_dict[remove_decoy_identification_str] = remove_decoy_identification
    for key, value in all_config_dict.items():
        if key == (pep_iden_str + fasta_database_str):
            config_dict[fasta_database_str] = value
        elif key == (pro_iden_str + decoy_prefix_str):
            config_dict[decoy_prefix_str] = value
        elif key == (pro_iden_str + min_peptide_per_protein_str):
            config_dict[min_peptide_per_protein_str] = value
        elif key == (pro_iden_str + min_unique_peptide_per_protein_str):
            config_dict[min_unique_peptide_per_protein_str] = value
        elif key == (pro_iden_str + remove_decoy_identification_str):
            config_dict[remove_decoy_identification_str] = value
        elif key == (pep_iden_str + cleave_after_residues_str):
            config_dict[cleave_after_residues_str] = value
        elif key == (pep_iden_str + cleave_before_residues_str):
            config_dict[cleave_before_residues_str] = value
        elif key == (pro_iden_str + FDR_Filtering_str):
            config_dict[FDR_Filtering_str] = value
        else:
            continue

    # return config dictionary
    return config_dict


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


def generate_psm_pep_txt(input_file, out_folder, psm_filtered_list):
    # get the FDR, # target, # decoy
    psm_target_int = 0
    psm_decoy_int = 0
    pep_target_int = 0
    pep_decoy_int = 0
    pep_set = set()
    for oPsm in psm_filtered_list:
        if not (oPsm.RealLabel == LabelFwd or oPsm.RealLabel == LabelTest):
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
    base_out_filename = input_file.split('/')[-1]
    base_out = out_folder + base_out_filename
    with open(base_out + ".psm.txt", 'w') as fw:
        
        fw.write('#    ########################################\n')
        fw.write('#    ####### PSM Filtering by Sipros ########\n')
        fw.write('#    ########################################\n')
        fw.write('#    \n')
        fw.write('#    #######################\n')
        fw.write('#    # Statistical Results #\n')
        fw.write('#    #######################\n')
        fw.write('#    \n')
        fw.write('#    [Statistical_Results]\n')
        fw.write('#    \n')
        fw.write('#    # Numbers of psm after filtering\n')
        fw.write('#    Decoy_PSMs_After_Filtering = %d\n' % psm_decoy_int)
        fw.write('#    Target_PSMs_After_Filtering = %d\n' % psm_target_int)
        fw.write('#    # PSM FDR = Decoy_PSMs_After_Filtering / Target_PSMs_After_Filtering\n')
        fw.write('#    PSM_FDR = %.2f%%\n' % (100.0*psm_decoy_int/psm_target_int))
        fw.write('#    \n')
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
            fw.write('Sipros10')
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
        if not (oPsm.RealLabel == LabelFwd or oPsm.RealLabel == LabelTest):
            continue
        pep_ID = oPsm.IdentifiedPeptide + '_+_' + str(oPsm.ParentCharge)
        if pep_ID in pep_sub_dict:
            pep_sub_dict[pep_ID].add(oPsm)
        else:
            oPeptide = Peptide()
            oPeptide.set(oPsm)
            pep_sub_dict[pep_ID] = oPeptide
    
    with open(base_out + ".pep.txt", 'w') as fw:
        # statistic results
        fw.write('#    ########################################\n')
        fw.write('#    ####### PSM Filtering by Sipros ########\n')
        fw.write('#    ########################################\n')
        fw.write('#    \n')
        fw.write('#    #######################\n')
        fw.write('#    # Statistical Results #\n')
        fw.write('#    #######################\n')
        fw.write('#    \n')
        fw.write('#    [Statistical_Results]\n')
        fw.write('#    \n')
        fw.write('#    # Numbers of peptide after filtering\n')
        fw.write('#    Decoy_peptides_After_Filtering = %d\n' % pep_decoy_int)
        fw.write('#    Target_peptides_After_Filtering = %d\n' % pep_target_int)
        fw.write('#    # peptide FDR = Decoy_peptides_After_Filtering / Target_peptides_After_Filtering\n')
        fw.write('#    peptide_FDR = %.2f%%\n' % (100.0*pep_decoy_int/pep_target_int))
        fw.write('#    \n')
        
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


## check pep and psm files pair set, and save run#
def get_run_num(pep_file_list, psm_file_list):

    # dictionary of Run# for each pep_file
    run_num_dict = {}
    psm_run_num_dict = {}

    # read multiple pep files
    for pep_file_idx, pep_file in enumerate(pep_file_list):
    
        # check file exist and open
        check_file_exist(pep_file)

        # If base common prefix ends with '.pep.txt', then remove '.pep.txt'
        base_pep_file = pep_file.replace(pep_file_ext, "")

        # make a psm filename using pep filename
        psm_file = base_pep_file + psm_file_ext

        # check psm file exist
        check_file_exist(psm_file)

        # for Run#
        run_tag = 'Run' + str(pep_file_idx + 1)

        # run_num_dict
        run_num_dict[pep_file] = run_tag
        psm_run_num_dict[psm_file] = run_tag

    return (run_num_dict, psm_run_num_dict)

PepOutFields = sipros_post_module.PepOutFields
# # need pep.txt file and pro.txt file, and big PSM table file
# # Spectral Count for original/identified peptide
# # Unique peptide counts and total peptide counts for protein
def feature_update(run_num_dict, pro_file, psm_tab, psm_tab_new):
    
    # save the pep file and pro file data to the defaultdict
    id_pep_data_dict = {}
    or_pep_data_dict = {}
    pro_data_dict = {}
    
    # key = pep_file, val = run_num , sorted by Run# index
    for pep_file, _run_num in sorted(run_num_dict.items(), key=lambda x: x[1][-1]):
        f = open(pep_file, 'rb')
        # read line with csv
        pep_reader = csv.reader(CommentedFile(f),
                                   delimiter='\t')
        # skip header
        _headline = pep_reader.next()
        
        # get data
        for pep_line in pep_reader:
            # adapt class PepOutFields
            pep_obj = PepOutFields._make(pep_line)
            identified_peptide = pep_obj.IdentifiedPeptide
            identified_peptide = identified_peptide.strip()
            ParentCharge = pep_obj.ParentCharge
            identified_peptide = identified_peptide +'_' + ParentCharge
            original_peptide = pep_obj.OriginalPeptide
            original_peptide = original_peptide.strip()
            
            iSpectralCount = int(pep_obj.SpectralCount.strip())
            
            if identified_peptide in id_pep_data_dict:
                id_pep_data_dict[identified_peptide] += iSpectralCount
            else:
                id_pep_data_dict[identified_peptide] = iSpectralCount
            
            if original_peptide in or_pep_data_dict:
                or_pep_data_dict[original_peptide] += iSpectralCount
            else:
                or_pep_data_dict[original_peptide] = iSpectralCount
        f.close()
    
    f = open(pro_file, 'rb')
    pro_reader = csv.reader(CommentedFile(f),
                                   delimiter='\t')
    # skip header
    headline = pro_reader.next()
    iNumRuns = (len(headline) - 2)/6
    iUniquePeptideCounts = 0
    iTotalPeptideCounts = 0
    for asWords in pro_reader:
        identified_protain = asWords[0]
        iUniquePeptideCounts = 0
        iTotalPeptideCounts = 0
        for i in range(iNumRuns):
            iUniquePeptideCounts += int(asWords[i*6+1])
            iTotalPeptideCounts += int(asWords[i*6+2])
        pro_data_dict[identified_protain] = [iUniquePeptideCounts, iTotalPeptideCounts]
    
    f.close()
    
    fr = open(psm_tab, 'rb')
    psm_reader = csv.reader(CommentedFile(fr),
                                   delimiter='\t')
    
    # skip header
    headline = psm_reader.next()
    
    with open(psm_tab_new, 'w') as f:
        # header
        f.write('\t'.join(headline))
        f.write('\tpep_psm\t')
        f.write('pro_pep')
        f.write('\n')
        
        for asWords in psm_reader:
            original_peptide = asWords[7]
            if original_peptide in or_pep_data_dict:
                iNumPepPsm = or_pep_data_dict[original_peptide]
                iNumPepPsm -= 1
                if iNumPepPsm < 0:
                    print('Error')
                    exit(1)
            else:
                iNumPepPsm = 0
            sProteins = asWords[12]
            asProteins = (sProteins[1:-1]).split()
            iNumTotalPep = 0
            iNumUniqPep = 0
            for sProtein in asProteins:
                if sProtein in pro_data_dict:
                    l = pro_data_dict[sProtein]
                    iNumTotalPep += l[1]
                    iNumUniqPep += l[0]
            f.write('\t'.join(asWords))
            f.write('\t')
            f.write(str(iNumPepPsm))
            f.write('\t')
            if iNumUniqPep > 1:
                f.write('2')
            elif iNumTotalPep > 1:
                f.write('1')
            else:
                f.write('0')
            f.write('\n')
    fr.close()
    
def clean_folder(output_folder):
    for the_file in os.listdir(output_folder):
        file_path = os.path.join(output_folder, the_file)
        try:
            if os.path.isfile(file_path):
                os.unlink(file_path)
            # elif os.path.isdir(file_path): shutil.rmtree(file_path)
        except Exception as e:
            print(e)


def show_measured_predicted_rt(psm_list, filename_prefix):
    fw_fwr = open(filename_prefix+"_fwr.txt", 'w')
    fw_fwr.write("measuread\tpredicted\n")
    fw_shu = open(filename_prefix+"_shu.txt", 'w')
    fw_shu.write("measuread\tpredicted\n")
    
    for oPsm in psm_list:
        if oPsm.RealLabel == LabelFwd:
            fw_fwr.write(str(oPsm.fRtMeasured))
            fw_fwr.write('\t')
            fw_fwr.write(str(oPsm.fRtPredict))
            fw_fwr.write('\n')
        else:
            fw_shu.write(str(oPsm.fRtMeasured))
            fw_shu.write('\t')
            fw_shu.write(str(oPsm.fRtPredict))
            fw_shu.write('\n')
    
    fw_fwr.close()
    fw_shu.close()
    
def stacking(block_num_int, psm_list):
    # # construct training data
    psm_filtered_list = []
    data_list = []
    label_list = []
    unknown_list = []
    id_list = []     
    for oPsm in psm_list:
        if len(oPsm.feature_list) == 0:
            print('check')
        unknown_list.append(oPsm.feature_list)
        id_list.append(oPsm.iInnerId)
        if oPsm.TrainingLabel != LabelUnknown:
            data_list.append(oPsm.feature_list)
            label_list.append(oPsm.TrainingLabel)
                
    data_np = np.array(data_list)
    sys.stdout.write(str(len(data_list)) + "\n")
    label_np = np.array(label_list)
    
    # split the data into B blocks
    train_data_block = [[], []]
    train_label_block = [[], []]
    block_idx_list = np.random.randint(0, block_num_int, size=len(data_np))
    data_idx_int = 0
    for c in block_idx_list:
        if c == 0:
            train_data_block[0].append(data_np[data_idx_int])
            train_label_block[0].append(label_np[data_idx_int])
        else:
            train_data_block[1].append(data_np[data_idx_int])
            train_label_block[1].append(label_np[data_idx_int])
        data_idx_int += 1
        
    # features and models
    feature_matrix = [[1, 2, 3, 15, 16, 17],    # score section
                      [5, 29],                  # mass section
                      [24],                     # digestion section
                      [25, 26, 27, 28],         # pep pro support section
                      [30]]                     # PTM section
    all_feature_list = [1, 2, 3, 5, 15, 16, 17, 24, 25, 26, 27, 28, 29, 30]
    tier_1 = []
    for i in range(len(feature_matrix)):
        tier_1.append(linear_model.LogisticRegression(penalty='l2', 
                                             dual=False, 
                                             C=1, 
                                             fit_intercept=True, 
                                             intercept_scaling=1, 
                                             class_weight='balanced', 
                                             random_state=None, 
                                             solver='liblinear', 
                                             max_iter=100, 
                                             multi_class='ovr', 
                                             verbose=0, 
                                             warm_start=False, 
                                             n_jobs=-1))
        # tier_1[i].fit(np.array(train_data_block[1])[:, feature_matrix[i]], np.array(train_label_block[1]))
        tier_1[i].fit(data_np[:, feature_matrix[i]], label_np)
        predict_np = tier_1[i].predict(data_np[:, feature_matrix[i]])
        # show_TP_TN_FP_FN(label_np, predict_np)
    
    
    
    
    tier_2 = linear_model.LogisticRegression(penalty='l2', 
                                             dual=False, 
                                             C=1, 
                                             fit_intercept=True, 
                                             intercept_scaling=1, 
                                             class_weight='balanced', 
                                             random_state=None, 
                                             solver='liblinear', 
                                             max_iter=100, 
                                             multi_class='ovr', 
                                             verbose=0, 
                                             warm_start=False, 
                                             n_jobs=-1)
    # # train tier_2
    train_tier_2_np = []
    for i in range(len(feature_matrix)):
        predict_np = tier_1[i].predict_proba(data_np[:, feature_matrix[i]])
        train_tier_2_np.append(predict_np[:, 1])
    train_tier_2_np = np.transpose(np.array(train_tier_2_np))
    tier_2.fit(train_tier_2_np, label_np)
    '''
    
    # # re-train tier_1
    tier_1 = []
    for i in range(len(feature_matrix)):
        tier_1.append(linear_model.LogisticRegression(penalty='l2', 
                                             dual=False, 
                                             C=1, 
                                             fit_intercept=True, 
                                             intercept_scaling=1, 
                                             class_weight='balanced', 
                                             random_state=None, 
                                             solver='liblinear', 
                                             max_iter=100, 
                                             multi_class='ovr', 
                                             verbose=0, 
                                             warm_start=False, 
                                             n_jobs=-1))
        tier_1[i].fit(data_np[:, feature_matrix[i]], label_np)
    '''
    # # testing
    unknown_np = np.array(unknown_list)
    test_tier_2_np = []
    for i in range(len(feature_matrix)):
        predict_np = tier_1[i].predict_proba(unknown_np[:, feature_matrix[i]])
        test_tier_2_np.append(predict_np[:, 1])
    test_tier_2_np = np.transpose(np.array(test_tier_2_np))
    predict_np = tier_2.predict_proba(test_tier_2_np)
    
    #predict_np = logreg.predict_proba(train_data_np)
    for i in range(len(predict_np)):
        psm_list[id_list[i]].fPredictProbability = predict_np[i, 1]
    
    
    psm_new_list = re_rank(psm_list)
    del psm_list[:]
    psm_list.extend(psm_new_list)
    psm_filtered_list_local = show_Fdr(psm_list, None, fdr=0.01)
    # show_Fdr_charge(psm_list)
    psm_filtered_list.extend(psm_filtered_list_local)
    sys.stdout.write('\n')
        
    return psm_filtered_list

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

def mass_filter(psm_list):
    psm_new_list = []
    for oPsm in psm_list:
        if abs(oPsm.fMassDiff) > mass_tolerance:
            continue
        psm_new_list.append(oPsm)
    return psm_new_list

def main(argv=None):
    if argv is None:
        argv = sys.argv

    # parse options
    (input_file, config_file, output_folder, negative_file) = parse_options(argv)
    
    # get the configuration parameters
    config_dict = parse_config(config_file)

    # get the categorized data
    (psm_dict, psm_neg_list) = categorize(input_file, negative_file)
    psm_list = []
    for key, lPsm in psm_dict.items():
        psm_list.extend(lPsm)
    
    i = 0
    for oPsm in psm_list:
        oPsm.iInnerId = i
        i += 1
        
    remark_concensus(psm_list)
    psm_list = mass_filter(psm_list)
    
    # generate features
    # generate_Prophet_features(psm_list, config_dict)
    generate_Prophet_features_test(psm_list, config_dict)
    # generate_pep_pro_features(psm_list, config_dict)
    generate_Prophet_features_group(psm_list, config_dict)
    
    # q product filtering
    count_list = [0, 0, 0]
    for key, lPsm in psm_dict.items():
        final_cutoff_score = 0.0
        mark_training_label(lPsm, final_cutoff_score, count_list, sKey = None)

    print('Before Machine Learning:\n\t# FWD\t# REV\t# SHU')
    print('\t%d\t%d\t%d' % (count_list[0], count_list[1], count_list[2]))
    
    # set feature all PSMs
    for oPsm in psm_list:
        oPsm.get_feature_list()
        
    '''
    # ensamble learning
    psm_filtered_list = test_stacking_MOE(psm_list)
    # stacking(2, psm_list)
    print 'Done.'
    return
    '''
    # machine learning
    #test_random_forest(psm_dict, psm_list)
    del feature_selection_list[:]
    feature_selection_list.extend([1, 2, 3, 5, 15, 16, 17, 24, 26, 28]) #
    # feature_selection_list.extend([1, 2, 3, 5, 24, 25, 26, 27, 28, 29]) #
    # psm_filtered_list = logistic_regression(psm_dict, psm_list)
    # deep learning
    # psm_filtered_list = test_DeepLearning(psm_dict, psm_list, psm_neg_list, 0.01/3, None)
    
    psm_filtered_list = logistic_regression_no_category(psm_dict, psm_list, psm_neg_list, 0.005/3, None, output_folder, input_file, config_dict)
    # logistic_regression_2_LR(psm_dict, psm_list, psm_neg_list, 0.01, None, output_folder, input_file)
    # train bias test
    # psm_filtered_list = test_train_bias(psm_dict, psm_list, psm_neg_list, 0.01, None)
    # neural network
    # psm_filtered_list = svm_one_class_test(psm_dict, psm_list, psm_neg_list, 0.01, None)
    # random forest
    # psm_filtered_list = test_random_forest(psm_dict, psm_list, psm_neg_list, 0.01/3, None, output_folder, input_file)
    # naive bayes
    # psm_filtered_list = test_naive_bayes(psm_dict, psm_list, psm_neg_list, 0.01, None)
    # SVM
    # psm_filtered_list = test_SVM(psm_dict, psm_list, psm_neg_list, 0.01, None)
    # ada boost
    # psm_filtered_list = test_AdaBoost(psm_dict, psm_list, psm_neg_list, 0.01/3, None)
    # ensemble voting
    # psm_filtered_list = test_voting(psm_list, 0.01)
    # ensemble LR
    # psm_filtered_list = test_LR_ensamble(psm_dict, psm_list, psm_neg_list, 0.01/3, None, output_folder, input_file)
    print('Done.')
    #return
    generate_psm_pep_txt(input_file, output_folder, psm_filtered_list)
    # ensemble MoE
    # psm_filtered_list = test_MoE_ensamble(psm_list)
    
    print('Done.')
    # return
    # RT data
    # get_RT_data(psm_filtered_list, output_folder, psm_list)
    # get_RT_elude_data(psm_filtered_list, output_folder, psm_list)
    # (psm_selected_list, psm_left_list) = filter_Fdr(psm_list, 0.01)
    # (psm_selected_list, psm_left_list) = filter_Fdr2(psm_list, 0.01)
    # show_measured_predicted_rt(psm_selected_list, '/media/xgo/Seagate/Proteomics/Experiments/Angelo/Sipros10/RTime/check/fdrless0.01')
    # show_measured_predicted_rt(psm_left_list, '/media/xgo/Seagate/Proteomics/Experiments/Angelo/Sipros10/RTime/check/fdrlarger0.01')
    # psm_filtered_list = logistic_regression_no_category_rt(psm_dict, psm_list, psm_neg_list, 0.01, psm_left_list)
    
    return
    

if __name__ == '__main__':
    sys.exit(main())
