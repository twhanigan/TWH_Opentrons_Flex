import datetime
import subprocess
import pandas as pd
import pickle
from io import BytesIO
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

num_samples = 10
target_concentration = 1
final_volume = 0.5
# Construct today's date in the expected format
today_date = datetime.date.today().strftime("%y%m%d")
file_name = f"TWH_Plate_{today_date}.xlsx"
file_url = f"http://172.23.226.47:48888/files/TWH/{file_name}"

def fetch_data():
    try:
        # Fetch the file using curl via subprocess
        result = subprocess.run(["curl", "-s", file_url], capture_output=True, check=True)
        data = result.stdout

        # Read the Excel file content into pandas DataFrame
        df = pd.read_excel(BytesIO(data), skiprows=5, usecols="B:N",index_col=0)
        # Return the filename and serialized DataFrame
        return file_name, df
    except subprocess.CalledProcessError as e:
        return None, f"Failed to fetch the file: {e}"
    except Exception as ex:
        return None, f"Error reading the Excel file: {ex}"

file_name, df = fetch_data()
if file_name:
    # Print the filename
    print(file_name)
    # Print the serialized DataFrame (binary data will be printed)
print(df)

# Create a list of well names (A1 to H12)
well_names = [f"{row}{col}" for col in range(1, 13) for row in "ABCDEFGH"]

# Flatten the absorbance values into a single list
absorbance_values = df.values.flatten()
print(absorbance_values)

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

#write the output to the instrument jupyter notebook directory
filename = f"Protocol_output_{today_date}.csv"
directory = Path("./")
output_file_destination_path = directory.joinpath(filename)
normalized_samples.to_csv(output_file_destination_path)

# Plot the data and linear regression
plt.figure(figsize=(8, 6))
plt.scatter(samples_1_to_8['Protein Concentration (mg/mL)'], samples_1_to_8['Mean Absorbance'], color='blue', label='Data')
plt.plot(samples_1_to_8['Protein Concentration (mg/mL)'], y_pred, color='red', label=f'Linear regression (R² = {r_squared:.2f})')
plt.xlabel('Protein Concentration (mg/mL)')
plt.ylabel('Mean Absorbance')
plt.title('Linear Regression of Protein Concentration vs Absorbance')
plt.legend()

# Save the plot as an image
plot_image_path = directory.joinpath(f"linear_regression_plot_{today_date}.png")
plt.savefig(plot_image_path)

# Optionally, display the plot
plt.show()