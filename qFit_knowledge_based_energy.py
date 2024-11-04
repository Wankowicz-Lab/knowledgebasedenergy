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
    # Initialize an empty DataFrame to store results
    final_results = pd.DataFrame()

    # Loop through each row in chi_data
    for x in range(len(chi_data)):
        row = chi_data.iloc[[x]].copy()
        eval_angle = row['angle'].values[0]  # Get the angle value

        # Filter chi_bins for matching bin range
        matched_row = chi_bins[(chi_bins['bin min'] < eval_angle) & (chi_bins['bin max'] > eval_angle)]

        # If there's a matching bin, add the E value to the row
        if not matched_row.empty:
            row.loc[:, 'E'] = matched_row['E'].values[0]  # Use .loc for assignment
            final_results = pd.concat([final_results, row], ignore_index=True)

    return final_results

# Function to generate scatter plots for each unique AA_CHI combination
def plot_data(final_sorted_data):
    
    sns.set(style='whitegrid', palette='muted')
    
    # Define colors for qFit and non-qFit data points
    qfit_color = sns.color_palette("muted")[0]  # Blueish tone
    non_qfit_color = sns.color_palette("muted")[3]  # Reddish tone

    # Get unique amino acid and chi-angle combinations
    unique_aa_chi = final_sorted_data['AA_CHI'].unique()

    # Iterate over each unique AA_CHI to generate plots
    for aa_chi in unique_aa_chi:
        specific_aa_chi = final_sorted_data[final_sorted_data['AA_CHI'] == aa_chi]
        amino_acid, chi_angle = aa_chi.split('_')
        chi_angle = chi_angle.replace('CHI', 'χ')
        chi_angle_label = f"{amino_acid} {chi_angle}(°)" 

        # Separate data based on the source file
        qFit = specific_aa_chi[specific_aa_chi['source_file'] == 'qFit_rotamers_output']
        non_qFit = specific_aa_chi[specific_aa_chi['source_file'] == 'rotamers_output']
        
        # Create the plot
        plt.figure(figsize=(10, 6))
        
        sns.scatterplot(
            data=non_qFit, 
            x='angle', 
            y='E', 
            color=non_qfit_color,  
            s=100, 
            edgecolor='gray', 
            alpha=0.9,
            label='Rotamers Output'
        )
        sns.scatterplot(
            data=qFit, 
            x='angle', 
            y='E', 
            color=qfit_color,  
            s=100, 
            edgecolor='gray', 
            alpha=0.9,
            label='qFit Rotamers Output'
        )
        
        plt.title(f'Knowledge-based ΔE ({amino_acid} {chi_angle})', fontsize=16, fontweight='bold')
        plt.xlabel(chi_angle_label, fontsize=14, fontweight='bold')
        plt.ylabel('ΔE (kcal/mol)', fontsize=14, fontweight='bold')
        plt.xlim(0, 360)
        plt.xticks([0, 120, 240, 360], fontsize=12)
        plt.yticks(fontsize=12)
        plt.grid(True, linestyle='--', alpha=0.5)
        plt.legend(title='Source File', title_fontsize='13', fontsize='12', loc='upper right')
        plt.tight_layout()
        plt.savefig(f'qFit_knowledge_based_energy{amino_acid}_{chi_angle}.png')

def main():
    qFit_b_factor_files = [os.path.join(qFit_rotamers, item) for item in os.listdir(qFit_rotamers) if item.endswith('qFit_rotamers_output.csv')]

    # Load DataFrames
    qFit_b_factor_dfs = []
    qFit_not_parsed = []

    for file in qFit_b_factor_files:
        try:
            qFit_b_factor_dfs.append(pd.read_csv(file))
        except:
            qFit_not_parsed.append(file)

    # Check and concatenate DataFrames if any were read
    if qFit_b_factor_dfs:
        combined_qFit_b_factor_df = pd.concat(qFit_b_factor_dfs, ignore_index=True)
        
        # Rename columns
        combined_qFit_b_factor_df.rename(columns={'residue_name': 'residue_type', 'rotamer_value': 'angle', 'nchi': 'chi_angle'}, inplace=True)
        
        # Ensure 'angle' column is numeric, setting errors='coerce' to handle non-numeric values by setting them to NaN
        combined_qFit_b_factor_df['angle'] = pd.to_numeric(combined_qFit_b_factor_df['angle'], errors='coerce')
        
        # Make all angles positive (only on non-NaN values)
        combined_qFit_b_factor_df['angle'] = combined_qFit_b_factor_df['angle'].where(combined_qFit_b_factor_df['angle'].isna(), combined_qFit_b_factor_df['angle'] % 360)
        
        # Update 'chi_angle' values
        combined_qFit_b_factor_df['chi_angle'] = combined_qFit_b_factor_df['chi_angle'].map({
            0: 'chi1',
            1: 'chi2',
            2: 'chi3',
            3: 'chi4'
        }).fillna(combined_qFit_b_factor_df['chi_angle'])

        # Handle 'altloc' column if it exists
        if 'altloc' in combined_qFit_b_factor_df.columns:
            combined_qFit_b_factor_df['residue_altloc'] = combined_qFit_b_factor_df['residue'].astype(str) + '_' + combined_qFit_b_factor_df['altloc'].fillna('')
            subset_with_residue_altloc = combined_qFit_b_factor_df[combined_qFit_b_factor_df['altloc'].notna()]
        else:
            subset_with_residue_altloc = pd.DataFrame()  # Empty DataFrame if 'altloc' column not found
        
    b_factor_files = [os.path.join(rotamers, item) for item in os.listdir(rotamers) if item.endswith('rotamers_output.csv')]

    # Load DataFrames
    b_factor_dfs = []
    not_parsed = []

    for file in b_factor_files:
        try:
            b_factor_dfs.append(pd.read_csv(file))  # Append DataFrame to b_factor_dfs, not b_factor_files
        except Exception as e:
            print(f"Could not parse {file}: {e}")
            not_parsed.append(file)

    # Check and concatenate DataFrames if any were read
    if b_factor_dfs:
        combined_b_factor_df = pd.concat(b_factor_dfs, ignore_index=True)

        # Rename columns
        combined_b_factor_df.rename(columns={'residue_name': 'residue_type', 'rotamer_value': 'angle', 'nchi': 'chi_angle'}, inplace=True)
        
        # Ensure 'angle' column is numeric, setting errors='coerce' to handle non-numeric values by setting them to NaN
        combined_b_factor_df['angle'] = pd.to_numeric(combined_b_factor_df['angle'], errors='coerce')
        
        # Make all angles positive (only on non-NaN values)
        combined_b_factor_df['angle'] = combined_b_factor_df['angle'].where(combined_b_factor_df['angle'].isna(), combined_b_factor_df['angle'] % 360)
        
        # Update 'chi_angle' values
        combined_b_factor_df['chi_angle'] = combined_b_factor_df['chi_angle'].map({
            0: 'chi1',
            1: 'chi2',
            2: 'chi3',
            3: 'chi4'
        }).fillna(combined_b_factor_df['chi_angle'])

    # filter and concatenation based on matching residues in qFit subset
    if not subset_with_residue_altloc.empty:
        # Get unique residues to match
        matching_residues = subset_with_residue_altloc['residue'].unique()
        
        # Filter combined DataFrame to retain only matching residues
        filtered_new_df = combined_b_factor_df[combined_b_factor_df['residue'].isin(matching_residues)]
        
        # Add source_file column in a single vectorized operation before concatenation
        subset_with_residue_altloc = subset_with_residue_altloc.assign(source_file='qFit_rotamers_output')
        filtered_new_df = filtered_new_df.assign(source_file='rotamers_output')
        
        # Concatenate the DataFrames
        concatenated_df = pd.concat([subset_with_residue_altloc, filtered_new_df], ignore_index=True)
    else:
        concatenated_df = pd.DataFrame()  # Empty DataFrame if no matching residues

    # List all CSV files in the energy bin folder
        csv_files = [file for file in os.listdir(folder_path) if file.endswith('.csv')]

        # Initialize an empty list to collect the DataFrames
        dfs = []

        # Loop through all CSV files and append the processed DataFrames to the list
        for file_name in csv_files:
            df = pd.read_csv(os.path.join(folder_path, file_name))
            df["AA_CHI"] = file_name.split('.')[0]  # Add the new column for AA_CHI & labeling each row by file name
            dfs.append(df)

        # Concatenate all DataFrames in the list into one DataFrame
        energy_bin = pd.concat(dfs, ignore_index=True)

        # Create an empty list to store processed data
        process_list = []

        # Iterate through each unique amino acid and chi angle, and process data accordingly
        for amino_acid in concatenated_df['residue_type'].unique():
            # Filter the DataFrame by the current amino acid
            temp_df = cconcatenated_df[concatenated_df['residue_type'] == amino_acid]
            
            for chi_angle in temp_df['chi_angle'].unique():
                # Ensure chi_angle is a string and contains 'chi'
                if isinstance(chi_angle, str) and 'chi' in chi_angle:
                    # Capitalize chi_angle if valid
                    chi_angle_capital = chi_angle.upper()
                    
                    # Create the energy_bins key
                    energy_bins = f"{amino_acid}_{chi_angle_capital}"
                    
                    # Filter the energy_bin DataFrame by the energy_bins key
                    temp_energy = energy_bin[energy_bin['AA_CHI'] == energy_bins]
                    
                    # Further filter temp_df by the current chi_angle
                    temp_df2 = temp_df[temp_df['chi_angle'] == chi_angle]
                    
                    # Process the chi data and append to process_list
                    process = process_chi_data(temp_df2, temp_energy)

                    # Add the AA_CHI column using the energy_bins key
                    process["AA_CHI"] = energy_bins
                    
                    process_list.append(process)

        # Concatenate all processed DataFrames into one final DataFrame
        final_sorted_data = pd.concat(process_list, ignore_index=True)

        # Plotting
        plot_data(final_sorted_data)

# Run the main function
if __name__ == '__main__':
    main()    
