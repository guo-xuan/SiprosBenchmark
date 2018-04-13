'''
Created on Nov 5, 2017

@author: naux
'''

import os

def spe2pep_to_sip(spe2pep_file_s, sip_file_s):
    with open(sip_file_s, 'w') as fw:
        with open(spe2pep_file_s, 'r') as fr:
            for line_s in fr:
                if line_s.startswith('#'):
                    fw.write(line_s)
                if line_s.startswith('+'):
                    spe_l = line_s.split()
                    i = 1
                if line_s.startswith('*'):
                    pep_l = line_s.split()
                    fw.write(spe_l[1]) #Filename
                    fw.write('\t')
                    fw.write(spe_l[2]) #ScanNumber
                    fw.write('\t')
                    fw.write(spe_l[3]) #ParentCharge
                    fw.write('\t')
                    fw.write(spe_l[4]) #MeasuredParentMass
                    fw.write('\t')
                    fw.write(pep_l[3]) #CalculatedParentMass
                    fw.write('\t')
                    fw.write(spe_l[5]) #ScanType
                    fw.write('\t')
                    fw.write(spe_l[6]) #SearchName
                    fw.write('\t')
                    fw.write('WeightSum') #ScoringFunction
                    fw.write('\t')
                    fw.write(str(i)) #Rank
                    i += 1
                    fw.write('\t')
                    fw.write(pep_l[4]) #Score
                    fw.write('\t')
                    fw.write(pep_l[1]) #IdentifiedPeptide
                    fw.write('\t')
                    fw.write(pep_l[2]) #OriginalPeptide
                    fw.write('\t')
                    fw.write(pep_l[-1]) #ProteinNames
                    fw.write('\n')
                          
def spe2pep_to_sip_folder(input_folder, output_folder):
    for f in os.listdir(input_folder):
        if f.endswith(".Spe2Pep.txt"):
            spe2pep_to_sip(os.path.join(input_folder, f), os.path.join(output_folder, f).replace('Spe2Pep.txt', 'sip'))

if __name__ == '__main__':
    input_folder = "/home/naux/temp/"
    output_folder = "/home/naux/temp/tmp/"
    spe2pep_to_sip_folder(input_folder, output_folder)
    print('Done.')
    