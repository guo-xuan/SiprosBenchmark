'''
Created on Jul 6, 2017

@author: xgo
'''
import sys
from os import listdir
from os.path import isfile, join

def get_id_dict(original_folder_str, renamed_folder_str):
    file_list = [f for f in listdir(renamed_folder_str) if isfile(join(renamed_folder_str, f))]
    rtime_mass_scannumber_dict = {}
    scannumber_str = ''
    mass_str = ''
    time_str = ''

    for one_file in file_list:
        with open(join(renamed_folder_str, one_file), 'r') as fr:
            for line_str in fr:
                if line_str.startswith('S'):
                    split_list = line_str.strip().split('\t')
                    scannumber_str = split_list[1]
                    mass_str = split_list[3]
                if line_str.startswith('I\tRTime'):
                    split_list = line_str.strip().split('\t')
                    time_str = split_list[2]
                    unique_id_str = time_str + '_' + mass_str
                    if unique_id_str not in rtime_mass_scannumber_dict:
                        rtime_mass_scannumber_dict[unique_id_str] = scannumber_str
                    else:
                        print('error')
    
    id_dict = {}
    file_list = [f for f in listdir(original_folder_str) if isfile(join(original_folder_str, f))]
    for one_file in file_list:
        file_id_str = one_file[one_file.rfind('Run')+3:one_file.rfind('.')]
        with open(join(original_folder_str, one_file), 'r') as fr:
            for line_str in fr:
                if line_str.startswith('S'):
                    split_list = line_str.strip().split('\t')
                    scannumber_str = split_list[1]
                    mass_str = split_list[3]
                if line_str.startswith('I\tRTime'):
                    split_list = line_str.strip().split('\t')
                    time_str = split_list[2]
                    unique_id_str = time_str + '_' + mass_str
                    if unique_id_str not in rtime_mass_scannumber_dict:
                        print('something wrong.')
                        print(one_file)
                        print(unique_id_str)
                        sys.exit(1)
                    else:
                        id_dict[file_id_str + '_' + str(scannumber_str)] = rtime_mass_scannumber_dict[unique_id_str]
    
    print('get_id_dict is done.')
    return id_dict

def get_sipros_psms(id_dict, input_folder_str, output_file_str):
    fw = open(output_file_str, 'w')
    
    file_list = [f for f in listdir(input_folder_str) if f.endswith('psm.txt') and isfile(join(input_folder_str, f))]
    for one_file in file_list:
        with open(join(input_folder_str, one_file), 'r') as fr:
            for line_str in fr:
                if line_str.startswith('#') or line_str.startswith('Filename'):
                    continue
                split_list = line_str.strip().split('\t')
                if split_list[-1] == 'F':
                    continue
                file_id_str = split_list[0][split_list[0].rfind('Run')+3:split_list[0].rfind('.')]
                id_str = file_id_str + '_' + split_list[1]
                if id_str not in id_dict:
                    print('not in the dictionary')
                    continue
                new_id_str = id_dict[id_str]
                fw.write(new_id_str + '_' + split_list[2])
                fw.write('\t')
                fw.write(split_list[14])
                fw.write('\n')
    fw.close()
    print('get_sipros_psms is done.')
    
def get_psm_from_comet(input_file_str, output_file_str):
    with open(output_file_str, 'w') as fw:
        with open(input_file_str, 'r') as fr:
            fr.readline()
            for line_str in fr:
                split_list = line_str.strip().split('\t')
                if split_list[-1] == 'F':
                    continue
                scannumber_str = split_list[1]
                charge_str = split_list[2]
                pep_str = split_list[14]
                fw.write(scannumber_str+'_'+charge_str)
                fw.write('\t')
                fw.write(pep_str)
                fw.write('\n')
    
    print('get_psm_from_comet is done.')
    
def get_venn_number():
    file_list = ['/media/xgo/Seagate/Proteomics/Experiments/BenchmarkRev/Sipros10/0820/venn_psm/siproensemble.txt',
                 '/media/xgo/Seagate/Proteomics/Experiments/BenchmarkRev/Sipros10/0820/venn_psm/comet.txt',
                 '/media/xgo/Seagate/Proteomics/Experiments/BenchmarkRev/Sipros10/0820/venn_psm/myrimatch.txt',
                 '/media/xgo/Seagate/Proteomics/Experiments/BenchmarkRev/Sipros10/0820/venn_psm/msgf.txt']
    set_list = [set(), set(), set(), set()]
    for idx, one_file in enumerate(file_list):
        with open(one_file, 'r') as fr:
            for line_str in fr:
                line_str = line_str.strip()
                set_list[idx].add(line_str)
    
    return set_list
    
def get_overlap_number(set_list, set_selected_list):
    first_set = set_selected_list[0]
    rest_list = set_selected_list[1:]
    rest_set_list = []
    for i in rest_list:
        rest_set_list.append(set_list[i])
    count_int = 1
    total_int = len(set_selected_list)
    overlap_number_int = 0
    for psm_str in set_list[first_set]:
        count_int = 1
        for one_set in rest_set_list:
            if psm_str in one_set:
                count_int += 1
        
        if count_int == total_int:
            overlap_number_int+= 1
    
    print(str(overlap_number_int))
    return overlap_number_int
    
def draw_venn_diagram():
    pass
                        
def main():
    '''
    # extract psm from sipros ensemble
    original_folder_str = '/media/xgo/Seagate/Proteomics/Data/Benchmark/08202012/ms2/'
    renamed_folder_str = '/media/xgo/Seagate/Proteomics/Data/Benchmark/08202012/ms2_renamed_split/'
    id_dict = get_id_dict(original_folder_str, renamed_folder_str)
    
    input_folder_str = '/media/xgo/Seagate/Proteomics/Experiments/BenchmarkRev/Sipros10/0820/fdr_0.01/'
    output_file_str = '/media/xgo/Seagate/Proteomics/Experiments/BenchmarkRev/Sipros10/0820/venn_psm/siproensemble.txt'
    get_sipros_psms(id_dict, input_folder_str, output_file_str)
    '''
    '''
    # extract psm from comet
    input_file_str = '/media/xgo/Seagate/Proteomics/Experiments/BenchmarkRev/Myrimatch_4_percolator/0820/fdr_0.01/all.psm.txt'
    output_file_str = '/media/xgo/Seagate/Proteomics/Experiments/BenchmarkRev/Sipros10/0820/venn_psm/myrimatch.txt'
    get_psm_from_comet(input_file_str, output_file_str)
    '''
    '''
    # get the venn number
    set_list = get_venn_number()
    for one_set in set_list:
        print(len(one_set))
    set_selected_list = [0, 1]
    print(set_selected_list)
    get_overlap_number(set_list, set_selected_list)
    set_selected_list = [0, 2]
    print(set_selected_list)
    get_overlap_number(set_list, set_selected_list)
    set_selected_list = [0, 3]
    print(set_selected_list)
    get_overlap_number(set_list, set_selected_list)
    set_selected_list = [1, 2]
    print(set_selected_list)
    get_overlap_number(set_list, set_selected_list)
    set_selected_list = [1, 3]
    print(set_selected_list)
    get_overlap_number(set_list, set_selected_list)
    set_selected_list = [2, 3]
    print(set_selected_list)
    get_overlap_number(set_list, set_selected_list)
    set_selected_list = [0, 1, 2]
    print(set_selected_list)
    get_overlap_number(set_list, set_selected_list)
    set_selected_list = [0, 1, 3]
    print(set_selected_list)
    get_overlap_number(set_list, set_selected_list)
    set_selected_list = [1, 2, 3]
    print(set_selected_list)
    get_overlap_number(set_list, set_selected_list)
    set_selected_list = [0, 1, 2, 3]
    print(set_selected_list)
    get_overlap_number(set_list, set_selected_list)
    '''
if __name__ == '__main__':
    sys.exit(main())