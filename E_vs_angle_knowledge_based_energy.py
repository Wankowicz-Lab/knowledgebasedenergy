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
    # Set up aesthetics
    sns.set(style='whitegrid', palette='muted')
    
    # Get unique AA_CHI values
    unique_aa_chi = subset_data['AA_CHI'].unique()
    
    # Loop through each unique AA_CHI to create a scatter plot of ΔE vs Δangle
    for aa_chi in unique_aa_chi:
        # Filter data for the specific amino acid and chi angle
        specific_aa_chi = subset_data[subset_data['AA_CHI'] == aa_chi]
    
        # Extract amino acid and chi angle, replace 'CHI' with Greek letter 'χ'
        amino_acid, chi_angle = aa_chi.split('_')
        chi_angle = chi_angle.replace('CHI', 'χ')
        title_label = f"{amino_acid} {chi_angle}"
    
        # Create the plot
        plt.figure(figsize=(6, 6))  # Adjust figure to be square
        sns.scatterplot(data=specific_aa_chi, x='ΔE', y='Δangle', color='dodgerblue', s=100, edgecolor='w', alpha=0.8)
        plt.title(f'ΔE vs Δangle for {title_label}', fontsize=15, fontweight='bold')
        plt.xlabel('ΔE (kcal/mol)', fontsize=14)
        plt.ylabel(f'{chi_angle} Δangle (°)', fontsize=14)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.tight_layout()
        plt.savefig(f'E_vs_angle_knowledge_based_energy{amino_acid}_{chi_angle}.png')

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
        temp_df = concatenated_df[concatenated_df['residue_type'] == amino_acid]
        
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
    
    # Separate the qFit and rotamers output data
    qfit_data = final_sorted_data[final_sorted_data['source_file'] == 'qFit_rotamers_output']
    rotamers_data = final_sorted_data[final_sorted_data['source_file'] == 'rotamers_output']
    
    # Merge qFit and rotamers data on matching residue, chain, chi_angle, and AA_CHI for comparison
    merged_data = pd.merge(
        qfit_data, rotamers_data,
        on=['chain', 'residue', 'residue_type', 'chi_angle', 'AA_CHI'],
        suffixes=('_qFit', '_rotamers')
    )
    
    # Calculate ΔE and Δangle
    merged_data['ΔE'] = merged_data['E_qFit'] - merged_data['E_rotamers']
    merged_data['Δangle'] = merged_data['angle_qFit'] - merged_data['angle_rotamers']
    
    # Select relevant columns for the new subsetted data
    subset_data = merged_data[['chain', 'residue', 'residue_type', 'chi_angle', 'AA_CHI', 'ΔE', 'Δangle']]
    
    # Plotting
    plot_data(final_sorted_data)

# Run the main function
if __name__ == '__main__':
    main()
