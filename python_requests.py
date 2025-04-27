import requests
import pandas as pd
from io import BytesIO

# URL of the file (replace with the actual file URL)
file_url = "http://172.23.226.47:48888/files/TWH/TWH_Plate_250427.xlsx"

# Send the GET request to fetch the file
response = requests.get(file_url)

# Check if the request was successful
if response.status_code == 200:
    # Load the Excel file into a pandas DataFrame directly from the response content
    df = pd.read_excel(BytesIO(response.content), skiprows=6)
    df = df.drop(df.columns[0], axis=1)
    # Display the first few rows of the data
    print(df)
else:
    print(f"Failed to download file. Status code: {response.status_code}")