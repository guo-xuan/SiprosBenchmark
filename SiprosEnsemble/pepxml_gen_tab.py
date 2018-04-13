'''
Created on Sep 30, 2017

@author: naux
'''

def select_top_line_index(line_matrix, column_idx):
    idx = -1
    min_val = -1000000000
    for i, v in enumerate(line_matrix):
        if float(v[column_idx]) > min_val:
            min_val = float(v[column_idx])
            idx = i
    return idx

def output_record(tmp_list, fw, diff_idx, score_idx):
    charge_str = tmp_list[2]
    is_decoy = 'T'
    protein_name_list = tmp_list[12][1:-1].split(',')
    for protein in protein_name_list:
        if not 'Rev' in protein:
            is_decoy = 'F'
            
    diff = tmp_list[diff_idx]
    score = tmp_list[score_idx]
    psm_str = id_str + "\t" + charge_str + "\t" + score + "\t" + diff + "\t" + is_decoy + "\n"
    fw.write(psm_str)

if __name__ == '__main__':
    file_name = '/media/naux/naux_data/proteomics/experiments/pepxml/ecoli_full/EColi_Try_HCD_DE10ppm_CS_1000_NCE30_180_Run.tab'
    
    mvh_tab_file = '/media/naux/naux_data/proteomics/experiments/pepxml/ecoli_full/mvh_table.txt'
    xcorr_tab_file = '/media/naux/naux_data/proteomics/experiments/pepxml/ecoli_full/xcorr_table.txt'
    wdp_tab_file = '/media/naux/naux_data/proteomics/experiments/pepxml/ecoli_full/wdp_table.txt'
    
    fw_mvh = open(mvh_tab_file, 'w')
    fw_xcorr = open(xcorr_tab_file, 'w')
    fw_wdp = open(wdp_tab_file, 'w')
    
    fw_mvh.write("ID\tcharge\tmvh_score\tmvh_diff\tis_decoy\n")
    fw_xcorr.write("ID\tcharge\txcorr_score\txcorr_diff\tis_decoy\n")
    fw_wdp.write("ID\tcharge\twdp_score\twdp_diff\tis_decoy\n")
    
    all_psm_dic = {}
    
    with open(file_name, 'r') as fr:
        fr.readline()
        for line in fr:
            tmp_list = line.split()
            file_name_str = tmp_list[0].split('.')[0].split('_')[-1]
            scan_number_str = tmp_list[1]
            id_str = scan_number_str + "_" + file_name_str
            if id_str in all_psm_dic:
                all_psm_dic[id_str].append(line)
            else:
                all_psm_dic[id_str] = [line]
    
    for id_str, psm in all_psm_dic.items():
        line_matrix = []
        for tmp in psm:
            line_matrix.append(tmp.split())
        
        mvh_idx = select_top_line_index(line_matrix, 9)
        xcorr_idx = select_top_line_index(line_matrix, 10)
        wdp_idx = select_top_line_index(line_matrix, 11)
        
        output_record(line_matrix[mvh_idx], fw_mvh, 13, 9)
        output_record(line_matrix[xcorr_idx], fw_xcorr, 14, 10)
        output_record(line_matrix[wdp_idx], fw_wdp, 15, 11)
            
    fw_mvh.close()
    fw_xcorr.close()
    fw_wdp.close()
    
    print("Done.")
            