import os
import pandas as pd
from tqdm import tqdm
import seaborn as sns
import matplotlib.pyplot as plt

# Define the directories containing the files
rotamer = '/Users/joshua/Desktop/Wankowicz Lab/rotamer_output'
folder_path = '/Users/joshua/Desktop/Wankowicz Lab/backbone_independent_energy'

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

def main():
    # Use glob to find all *_B_factors.csv files in the specified directory
    b_factor_files = [os.path.join(rotamer,item) for item in os.listdir(rotamer)]
    
    # Initialize a list to store the DataFrames and track unparsed files
    b_factor_dfs = []
    not_parsed = []

    # Loop through each file and read it into a DataFrame
    for file in b_factor_files:
        try:
            df = pd.read_csv(file)
            # Optionally, add a new column to identify the source file
            b_factor_dfs.append(df)
        except:
            not_parsed.append(file)

    # Check if any DataFrames were added
    print("Number of DataFrames read:", len(b_factor_dfs))

    # Concatenate all DataFrames into a single DataFrame
    if b_factor_dfs:  # Check if the list is not empty
        combined_b_factor_df = pd.concat(b_factor_dfs, ignore_index=True)

        # Rename the columns
        combined_b_factor_df.rename(columns={'residue_name': 'residue_type', 'rotamer_value': 'angle', 'nchi': 'chi_angle'}, inplace=True)
        combined_b_factor_df['chi_angle'] = combined_b_factor_df['chi_angle'].replace({
            0: 'chi1',
            1: 'chi2',
            2: 'chi3',
            3: 'chi4'
        })

        # Make all angles positive
        combined_b_factor_df.loc[combined_b_factor_df['angle'] < 0, 'angle'] += 360    

    print("Not parsed:",len(not_parsed))
    
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
    for amino_acid in tqdm(combined_b_factor_df['residue_type'].unique()):
        # Filter the DataFrame by the current amino acid
        temp_df = combined_b_factor_df[combined_b_factor_df['residue_type'] == amino_acid]
        
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
    print(final_sorted_data)

    # Plotting
    plot_data(final_sorted_data)

# Function to plot data
def plot_data(final_sorted_data):
    sns.set(style='whitegrid', palette='muted') 

    # Get unique AA_CHI values
    unique_aa_chi = final_sorted_data['AA_CHI'].unique()

    # Loop through each unique AA_CHI
    for aa_chi in unique_aa_chi:
        # Filter for the specific amino acid and chi angle
        specific_aa_chi = final_sorted_data[final_sorted_data['AA_CHI'] == aa_chi]

        # Extract amino acid and chi angle, and replace 'CHI' with Greek letter 'χ'
        amino_acid, chi_angle = aa_chi.split('_')
        chi_angle = chi_angle.replace('CHI', 'χ')
        chi_angle_label = f"{amino_acid} {chi_angle}(°)" 

        # Create the plot
        plt.figure(figsize=(5, 3)) 
        sns.scatterplot(data=specific_aa_chi, x='angle', y='E', color='dodgerblue', s=100, edgecolor='w', alpha=0.8)
        plt.title(f'Knowledge-based ΔE({amino_acid} {chi_angle})', fontsize=15, fontweight='bold', color='gray')  # Set title color to gray
        plt.xlabel(chi_angle_label, fontsize=14)
        plt.ylabel('ΔE(kcal/mol)', fontsize=14)
        plt.xlim(0, 360)
        plt.xticks([0, 120, 240, 360], fontsize=12)
        plt.yticks(fontsize=12)
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.tight_layout()
        plt.show()

# Run the main function
if __name__ == '__main__':
    main()
