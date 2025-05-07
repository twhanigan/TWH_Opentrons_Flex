import pandas as pd
import numpy as np
#import matplotlib.pyplot as plt
import subprocess
from pathlib import Path
import datetime
import time

#Speed of pipetting NP40 lysis buffer=0.35, 2M Urea in EPPS=0.3
speed= 0.35
target_concentration = 1
final_volume = 0.440
final_volume_ul = final_volume*1000
num_samples = 10 #change this to the number of samples you need to run. The maximum is 18.
num_rows = 8  # A-H
num_replicates = 3  # the number of replicates
# Get today's date in YYMMDD format


# assign sample locations dynamically
sample_locations = []
for i in range(num_samples):
    if i < 6:  # B1 to B6
        sample_locations.append(f'B{i + 1}')
    elif i < 12:  # C1 to C6
        sample_locations.append(f'C{i - 5}')
    elif i < 18:  # D1 to D6
        sample_locations.append(f'D{i - 11}')
    else:
        break  # Stop if we exceed the number of available rows/columns


today_date = datetime.date.today().strftime("%y%m%d")

# Read the data file
df = pd.read_excel('WholeProt_3_250506.xlsx', header=5, nrows=8, usecols="C:N")

# Create a list of well names (A1 to H12)
well_names = [f"{row}{col}" for col in range(1, 13) for row in "ABCDEFGH"]

# Flatten the absorbance values into a single list
absorbance_values = df.values.flatten()

# Create the DataFrame
initial_df = pd.DataFrame({'Well': well_names, 'Absorbance': absorbance_values})

# Process data for normalization
samples, replicate_1, replicate_2, replicate_3 = [], [], [], []
sample_index = 1
for col_offset in range(0, 12, 3):  # Iterate by column groups (triplets)
    for row_offset in range(8):  # Iterate row-wise for each sample
        start = row_offset * 12 + col_offset  # Starting index for the sample
        if start + 2 < len(initial_df):
            samples.append(f"Sample {sample_index}")
            replicate_1.append(initial_df.iloc[start]['Absorbance'])
            replicate_2.append(initial_df.iloc[start + 1]['Absorbance'])
            replicate_3.append(initial_df.iloc[start + 2]['Absorbance'])
            sample_index += 1

final_df = pd.DataFrame({
    'Sample': samples,
    'Replicate 1': replicate_1,
    'Replicate 2': replicate_2,
    'Replicate 3': replicate_3
})

samples_1_to_8 = final_df.iloc[:8]
samples_1_to_8['Mean Absorbance'] = samples_1_to_8[['Replicate 1', 'Replicate 2', 'Replicate 3']].mean(axis=1)
protein_concentrations = [10, 5, 2.5, 1.25, 0.625, 0.3125, 0.15625, 0]
samples_1_to_8['Protein Concentration (mg/mL)'] = protein_concentrations

slope, intercept = np.polyfit(samples_1_to_8['Protein Concentration (mg/mL)'], samples_1_to_8['Mean Absorbance'], 1)
y_pred = slope * samples_1_to_8['Protein Concentration (mg/mL)'] + intercept
ss_res = np.sum((samples_1_to_8['Mean Absorbance'] - y_pred) ** 2)
ss_tot = np.sum((samples_1_to_8['Mean Absorbance'] - np.mean(samples_1_to_8['Mean Absorbance'])) ** 2)
r_squared = 1 - (ss_res / ss_tot)

unknown_samples = final_df.iloc[8:8 + num_samples]
unknown_samples['Mean Absorbance'] = unknown_samples[['Replicate 1', 'Replicate 2', 'Replicate 3']].mean(axis=1)
unknown_samples['Protein Concentration (mg/mL)'] = (unknown_samples['Mean Absorbance'] - intercept) / slope


unknown_samples['Sample Volume (mL)'] = (target_concentration * final_volume) / unknown_samples['Protein Concentration (mg/mL)']
unknown_samples['Diluent Volume (mL)'] = final_volume - unknown_samples['Sample Volume (mL)']
unknown_samples.loc[unknown_samples['Sample Volume (mL)'] > final_volume, ['Sample Volume (mL)', 'Diluent Volume (mL)']] = [final_volume, 0]
normalized_samples = unknown_samples[['Sample', 'Protein Concentration (mg/mL)', 'Sample Volume (mL)', 'Diluent Volume (mL)']].reset_index().drop(columns='index')

# Write the output and image of data plot to the instrument jupyter notebook directory
filename = f"Protocol_output_{today_date}.csv"
normalized_samples.to_csv(filename)

# Dilute sample in lysis buffer to 1 mg/ml on deep well plate
rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
destination_wells  = [f'{rows[i % 8]}{(i // 8)+ 1}' for i in range(len(normalized_samples))]
for i, row in normalized_samples.iterrows():
    source_well = sample_locations[i]
    normalized_volume = row['Sample Volume (mL)']
    diluent_volume = 500 - normalized_volume
    destination_well = destination_wells[i]
    #print(diluent_volume, 'reservoir[A7]', destination_well)
    #print(normalized_volume, source_well, destination_well)


# Transfer 10 ug in 10 uL volume to new column of plate3 and dilute with EPPS-Urea
tmt_dilution_wells = [f'{rows[i % 8]}{(i // 8) + 3}' for i in range(num_samples)]  # Assuming destination column is col 3
for i, row in normalized_samples.iterrows():
    protein_conc = row['Protein Concentration (mg/mL)']
    sample_vol_ul = min((10 / protein_conc), 10)
    diluent_vol = 10 - sample_vol_ul
    source = destination_wells[i]
    dest = tmt_dilution_wells[i]
    print(sample_vol_ul, source, dest)
    print(diluent_vol, 'reservoir[A12]', dest)