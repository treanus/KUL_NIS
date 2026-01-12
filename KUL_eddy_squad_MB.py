#!/usr/bin/env python
# our alternative version of eddy_squad
# Stefan Sunaert - 05/06/2024

import os
import sys
import json
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def read_json_files(folder_path):
    json_files = [pos_json for pos_json in os.listdir(folder_path) if pos_json.endswith('.json')]

    print(json_files)
    
    jsons_data = pd.DataFrame(columns=[
        'qc_mot_abs', 'qc_mot_rel', 'qc_cnr_avg', 'qc_cnr_std', 'qc_outliers_tot', 
        'qc_params_avg', 'qc_s2v_params_avg_std', 'qc_vox_displ_std', 'subject_id'
    ])

    for index, js in enumerate(json_files):
        with open(os.path.join(folder_path, js)) as json_file:
            json_text = json.load(json_file)
            qc_mot_abs = json_text.get('qc_mot_abs', None)
            qc_mot_rel = json_text.get('qc_mot_rel', None)
            qc_cnr_avg = json_text.get('qc_cnr_avg', None)
            qc_cnr_std = json_text.get('qc_cnr_std', None)
            qc_outliers_tot = json_text.get('qc_outliers_tot', None)
            qc_params_avg = json_text.get('qc_params_avg', None)
            qc_s2v_params_avg_std = json_text.get('qc_s2v_params_avg_std', None)
            qc_vox_displ_std = json_text.get('qc_vox_displ_std', None)
            subject_id = js.split('.')[0]  # Extract subject ID from the filename without extension

            jsons_data.loc[index] = [qc_mot_abs, qc_mot_rel, qc_cnr_avg, qc_cnr_std, qc_outliers_tot, 
                                     qc_params_avg, qc_s2v_params_avg_std, qc_vox_displ_std, subject_id]

    return jsons_data

def create_combined_violin_plot(data):
    plt.figure(figsize=(12, 10))
    print(data[['qc_outliers_tot']])

    # Plot for qc_mot_abs and qc_mot_rel
    plt.subplot(3, 1, 1)
    sns.violinplot(data=data[['qc_mot_abs', 'qc_mot_rel']], inner='quartiles')
    sns.stripplot(data=data[['qc_mot_abs', 'qc_mot_rel']], color='black', alpha=0.5, jitter=True, marker='o', size=5)
    plt.title('Distribution of qc_mot_abs and qc_mot_rel')
    plt.xlabel('Metrics')
    plt.ylabel('Values')
    plt.xticks(ticks=[0, 1], labels=['qc_mot_abs', 'qc_mot_rel'])
    
    
    # Plot for qc_cnr_avg
    qc_cnr_avg_list = list(map(list, zip(*data['qc_cnr_avg'].dropna())))
    print(list(map(list, zip(*qc_cnr_avg_list))))
    #for i, qc_cnr_avg in enumerate(qc_cnr_avg_list):
    for i, qc_cnr_avg in enumerate(qc_cnr_avg_list):
        print(i)
        plt.subplot(3, 4, i+5)
        sns.violinplot(data=[qc_cnr_avg], inner='quartiles')
        sns.stripplot(data=[qc_cnr_avg], color='black', alpha=0.5, jitter=True, marker='o', size=5)
        plt.title(f'Distribution of qc_cnr_avg {i+1}')
        plt.xlabel('qc_cnr_avg')
        plt.ylabel('Values')
    
    # Plot for qc_cnr_std
    qc_cnr_std_list = list(map(list, zip(*data['qc_cnr_std'].dropna())))
    for i, qc_cnr_std in enumerate(qc_cnr_std_list):
        print(i)
        plt.subplot(3, 4, i+5+len(qc_cnr_avg_list))
        sns.violinplot(data=[qc_cnr_std], inner='quartiles')
        sns.stripplot(data=[qc_cnr_std], color='black', alpha=0.5, jitter=True, marker='o', size=5)
        plt.title(f'Distribution of qc_cnr_std {i+1}')
        plt.xlabel('qc_cnr_std')
        plt.ylabel('Values')
    
    plt.tight_layout()
    plt.show()

def write_to_text_file(data, output_file):
    data.to_csv(output_file, index=False)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script_name.py <folder_path> <output_file.csv>")
        sys.exit(1)

    folder_path = sys.argv[1]
    output_file = sys.argv[2]

    json_data = read_json_files(folder_path)
    create_combined_violin_plot(json_data)
    write_to_text_file(json_data, output_file)
    print(f"Data written to {output_file}")
