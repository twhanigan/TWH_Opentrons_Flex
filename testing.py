import pandas as pd
import numpy as np
import subprocess
target_concentration = 1.5 #mg/ml
final_volume = 50 # uL
num_samples = 10 # change this to the number of samples you need to run. The maximum is 18.
num_rows = 8  # A-H
num_replicates = 3  # the number of replicates
speed= 0.35 #Speed of pipetting NP40 lysis buffer=0.35, 2M Urea in EPPS=0.3
df = pd.read_excel("C:\\Users\\thanigan\\Downloads\\OM_TWH_Exp5_250719(3).xlsx", header=5, nrows=8, usecols="C:N")

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
print(unknown_samples[['Sample', 'Protein Concentration (mg/mL)', 'Sample Volume (mL)', 'Diluent Volume (mL)']])

normalized_samples = unknown_samples[['Sample', 'Protein Concentration (mg/mL)', 'Sample Volume (mL)', 'Diluent Volume (mL)']].reset_index().drop(columns='index')
# Write the output and image of data plot to the instrument jupyter notebook directory

rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
destination_wells  = [f'{rows[i % 8]}{(i // 8)+ 1}' for i in range(len(normalized_samples))]



click_volume = 6*(final_volume/50)
print(click_volume)
loading_buffer_volume = round((final_volume+click_volume) / 3, 1)
print(loading_buffer_volume)
columns = sorted(set(well[1:] for well in destination_wells), key=int)
column_targets = [f'A{col}' for col in columns]
[print(well) for well in column_targets]

def mix_click_reagents():
    if volume_click_reaction < 100:
        positions_mixing = [1, 2, 3]
        pipette = p50_multi
    elif 100 < volume_click_reaction < 250:
        positions_mixing = [1, 4, 9]
        protocol.move_labware(labware=partial_50, new_location='B4', use_gripper=True)
        pipette = p1000_multi
    elif 250 < volume_click_reaction < 500:
        positions_mixing = [1, 6, 11]
        protocol.move_labware(labware=partial_50, new_location='B4', use_gripper=True)
        pipette = p1000_multi
    elif 500 < volume_click_reaction < 1000:
        positions_mixing = [1, 10, 16]
        protocol.move_labware(labware=partial_50, new_location='B4', use_gripper=True)
        pipette = p1000_multi
    else:
        positions_mixing = [1, 1, 1]  # Fallback/default case
        pipette = p1000_multi

    pipette.aspirate(final_volume/2, location_tube.bottom(z=positions_mixing[0]))
    pipette.dispense(final_volume/2, location_tube.bottom(z=positions_mixing[1]))
    pipette.aspirate(final_volume/3, location_tube.bottom(z=positions_mixing[2]))
    pipette.dispense(final_volume/3, location_tube.bottom(z=positions_mixing[0]))
    pipette.mix(7, volume_mixing, location_tube.bottom(z = position)) 
mix_click_reagents()
