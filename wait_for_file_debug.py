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

