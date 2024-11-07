import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Define the directories containing the files
rotamers = '/dors/wankowicz_lab/all_pdb/1_10000/output_rotamer/'
qFit_rotamers = '/dors/wankowicz_lab/all_pdb/1_10000/qFit_rotamers_output/'
folder_path = '/dors/wankowicz_lab/shared/backbone_independent_energy'

# Function to process chi data
def process_chi_data(chi_data, chi_bins): 
    final_results = pd.DataFrame()
    for x in range(len(chi_data)):
        row = chi_data.iloc[[x]].copy()
        eval_angle = row['angle'].values[0]
        matched_row = chi_bins[(chi_bins['bin min'] < eval_angle) & (chi_bins['bin max'] > eval_angle)]
        if not matched_row.empty:
            row.loc[:, 'E'] = matched_row['E'].values[0]
            final_results = pd.concat([final_results, row], ignore_index=True)
    return final_results

# Function to generate scatter plots
def plot_data(subset_data):
    sns.set(style='whitegrid', palette='muted')
    unique_aa_chi = subset_data['AA_CHI'].unique()
    for aa_chi in unique_aa_chi:
        specific_aa_chi = subset_data[subset_data['AA_CHI'] == aa_chi]
        amino_acid, chi_angle = aa_chi.split('_')
        chi_angle = chi_angle.replace('CHI', 'χ')
        title_label = f"{amino_acid} {chi_angle}"
        plt.figure(figsize=(6, 6))
        sns.scatterplot(data=specific_aa_chi, x='ΔE', y='Δangle', color='dodgerblue', s=100, edgecolor='w', alpha=0.8)
        plt.title(f'ΔE vs Δangle for {title_label}', fontsize=15, fontweight='bold')
        plt.xlabel('ΔE (kcal/mol)', fontsize=14)
        plt.ylabel(f'{chi_angle} Δangle (°)', fontsize=14)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.tight_layout()
        plt.savefig(f'E_vs_angle_knowledge_based_energy_{amino_acid}_{chi_angle}.png')

def main():
    qFit_b_factor_files = [os.path.join(qFit_rotamers, item) for item in os.listdir(qFit_rotamers) if item.endswith('qFit_rotamers_output.csv')]
    qFit_b_factor_dfs = []
    qFit_not_parsed = []

    for file in qFit_b_factor_files:
        try:
            qFit_b_factor_dfs.append(pd.read_csv(file))
        except:
            qFit_not_parsed.append(file)

    if qFit_b_factor_dfs:
        combined_qFit_b_factor_df = pd.concat(qFit_b_factor_dfs, ignore_index=True)
        combined_qFit_b_factor_df.rename(columns={'residue_name': 'residue_type', 'rotamer_value': 'angle', 'nchi': 'chi_angle'}, inplace=True)
        combined_qFit_b_factor_df['angle'] = pd.to_numeric(combined_qFit_b_factor_df['angle'], errors='coerce') % 360
        combined_qFit_b_factor_df['chi_angle'] = combined_qFit_b_factor_df['chi_angle'].map({0: 'chi1', 1: 'chi2', 2: 'chi3', 3: 'chi4'})

        if 'altloc' in combined_qFit_b_factor_df.columns:
            combined_qFit_b_factor_df['residue_altloc'] = combined_qFit_b_factor_df['residue'].astype(str) + '_' + combined_qFit_b_factor_df['altloc'].fillna('')
            subset_with_residue_altloc = combined_qFit_b_factor_df[combined_qFit_b_factor_df['altloc'].notna()]
        else:
            subset_with_residue_altloc = pd.DataFrame()

    b_factor_files = [os.path.join(rotamers, item) for item in os.listdir(rotamers) if item.endswith('rotamers_output.csv')]
    b_factor_dfs = []
    not_parsed = []

    for file in b_factor_files:
        try:
            b_factor_dfs.append(pd.read_csv(file))
        except Exception as e:
            print(f"Could not parse {file}: {e}")
            not_parsed.append(file)

    if b_factor_dfs:
        combined_b_factor_df = pd.concat(b_factor_dfs, ignore_index=True)
        combined_b_factor_df.rename(columns={'residue_name': 'residue_type', 'rotamer_value': 'angle', 'nchi': 'chi_angle'}, inplace=True)
        combined_b_factor_df['angle'] = pd.to_numeric(combined_b_factor_df['angle'], errors='coerce') % 360
        combined_b_factor_df['chi_angle'] = combined_b_factor_df['chi_angle'].map({0: 'chi1', 1: 'chi2', 2: 'chi3', 3: 'chi4'})

    if not subset_with_residue_altloc.empty:
        matching_residues = subset_with_residue_altloc['residue'].unique()
        filtered_new_df = combined_b_factor_df[combined_b_factor_df['residue'].isin(matching_residues)]
        subset_with_residue_altloc = subset_with_residue_altloc.assign(source_file='qFit_rotamers_output')
        filtered_new_df = filtered_new_df.assign(source_file='rotamers_output')
        concatenated_df = pd.concat([subset_with_residue_altloc, filtered_new_df], ignore_index=True)
    else:
        concatenated_df = pd.DataFrame()

    csv_files = [file for file in os.listdir(folder_path) if file.endswith('.csv')]
    dfs = []

    for file_name in csv_files:
        df = pd.read_csv(os.path.join(folder_path, file_name))
        df["AA_CHI"] = file_name.split('.')[0]
        dfs.append(df)

    energy_bin = pd.concat(dfs, ignore_index=True)
    process_list = []

    for amino_acid in concatenated_df['residue_type'].unique():
        temp_df = concatenated_df[concatenated_df['residue_type'] == amino_acid]
        
        for chi_angle in temp_df['chi_angle'].unique():
            if isinstance(chi_angle, str) and 'chi' in chi_angle:
                chi_angle_capital = chi_angle.upper()
                energy_bins = f"{amino_acid}_{chi_angle_capital}"
                temp_energy = energy_bin[energy_bin['AA_CHI'] == energy_bins]
                temp_df2 = temp_df[temp_df['chi_angle'] == chi_angle]
                process = process_chi_data(temp_df2, temp_energy)
                process["AA_CHI"] = energy_bins
                process_list.append(process)

    final_sorted_data = pd.concat(process_list, ignore_index=True)
    qfit_data = final_sorted_data[final_sorted_data['source_file'] == 'qFit_rotamers_output']
    rotamers_data = final_sorted_data[final_sorted_data['source_file'] == 'rotamers_output']
    
    merged_data = pd.merge(
        qfit_data, rotamers_data,
        on=['chain', 'residue', 'residue_type', 'chi_angle', 'AA_CHI'],
        suffixes=('_qFit', '_rotamers')
    )
    
    merged_data['ΔE'] = merged_data['E_rotamers'] - merged_data['E_qFit']
    merged_data['Δangle'] = merged_data['angle_rotamers'] - merged_data['angle_qFit']
    subset_data = merged_data[['chain', 'residue', 'residue_type', 'chi_angle', 'AA_CHI', 'ΔE', 'Δangle']]
    plot_data(subset_data)

if __name__ == '__main__':
    main()
