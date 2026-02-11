from opentrons import protocol_api
from opentrons.protocol_api import SINGLE, ALL
import pandas as pd
import numpy as np
import subprocess
from pathlib import Path
#import matplotlib.pyplot as plt
import datetime
import time
import re

metadata = {
    'protocolName': 'Normalization Western blotting',
    'author': 'Assistant',
    'description': 'BCA and normalization for western blotting'
}

requirements = {
    "robotType": "Flex",
    "apiLevel": "2.21"
}

def add_parameters(parameters):

    parameters.add_int(
        variable_name="num_samples",
        display_name="Number of samples",
        description="Number of input samples to be tested for mycoplasma.",
        default=8,
        minimum=1,
        maximum=24
    )
    parameters.add_int(
        variable_name="standards_col",
        display_name="standards column",
        description="Integer of column on plate1 where standards will be diluted",
        default=1,
        minimum=1,
        maximum=12
    )
    parameters.add_int(
        variable_name="ug_protein",
        display_name="µg of protein",
        description="Amount of lysate to load for western blot",
        default=40,
        minimum=5,
        maximum=100,
        unit="µg"
    )

    parameters.add_int(
        variable_name="final_volume",
        display_name="final volume lysate",
        description="Volume to normalize µg of protein in",
        default=20,
        minimum=10,
        maximum=100,
        unit="µL"
    )

def run(protocol: protocol_api.ProtocolContext):
    protocol.comment(
        "Place BSA Standard in A1, Lysis buffer in A2, tbta in A3, biotin in A4, cuso4 in A5, tcep in A6 and samples in row B")
    protocol.comment("Running the BCA assay")

    #Edit these
    target_concentration = protocol.params.ug_protein/protocol.params.final_volume # 40 ug protein/ 15 uL final volume
    speed= 0.35

    # Load modules
    heater_shaker = protocol.load_module('heaterShakerModuleV1', 'D1')
    thermocycler = protocol.load_module('thermocyclerModuleV2')
    temp_module = protocol.load_module('temperature module gen2', 'C1')
    mag_block = protocol.load_module('magneticBlockV1', 'D2')
    chute = protocol.load_waste_chute()

    # Load adapters
    temp_adapter = temp_module.load_labware('opentrons_24_aluminumblock_nest_1.5ml_snapcap')

    #set the heater_shaker temp to 60C
    heater_shaker.set_and_wait_for_temperature(50)

    #set the temp module to 0c
    temp_module.set_temperature(celsius=10)

    #open the thermocycler lid
    thermocycler.open_lid()
    
    # Load labware
    partial_50 = protocol.load_labware(load_name="opentrons_flex_96_filtertiprack_50ul",location="B3")
    tips_200 = protocol.load_labware(load_name="opentrons_flex_96_filtertiprack_200ul",location="A2")
    #tips_1000 = protocol.load_labware('opentrons_flex_96_filtertiprack_1000ul', 'C4')
    #plate1 = protocol.load_labware('opentrons_96_wellplate_200ul_pcr_full_skirt', 'A2')
    #plate2 = protocol.load_labware('corning_96_wellplate_360ul_flat', 'B2')
    plate3 = thermocycler.load_labware('nest_96_wellplate_100ul_pcr_full_skirt')    
    reservoir = protocol.load_labware('nest_12_reservoir_15ml', 'C2')
    
    # Liquid definitions
    bsa_standard = protocol.define_liquid(name = 'BSA Standard', display_color="#704848",)
    bsa_reag_a = protocol.define_liquid(name = 'Reagent A', display_color="#C0C0C0",)
    bsa_reag_b = protocol.define_liquid(name = 'Reagent B', display_color="#008000",)
    bsa_reag_c = protocol.define_liquid(name = 'Reagent C', display_color="#4B9CD3",)
    excess_lysis = protocol.define_liquid(name='Excess Lysis Buffer', display_color="#FFC0CB")  # Pink
    loading_buffer = protocol.define_liquid(name = 'Loading Buffer', display_color="#4169E1",)
    sample_liquids = [protocol.define_liquid(name = f'Sample {i + 1}', display_color="#FFA000",) for i in range(protocol.params.num_samples)]
    reduce_10X = protocol.define_liquid(name = 'Reducing Agent 10X', display_color="#4169E1",)
    empty_epp = protocol.define_liquid(name = 'Reducing Agent 10X', display_color="#FFFFFF",)

    # BSA and loading buffer
    temp_adapter['A1'].load_liquid(liquid=bsa_standard, volume=1000)  # 20 mg/ml BSA standard
    temp_adapter['A2'].load_liquid(liquid=loading_buffer, volume=1000)  # Additional lysis buffer for SP3
    temp_adapter['A3'].load_liquid(liquid=reduce_10X, volume=1000)  # Additional lysis buffer for SP3
    temp_adapter['A4'].load_liquid(liquid=empty_epp, volume=1000)  # Additional lysis buffer for SP3

    # Reservoir assignments for washes and digestion
    reservoir['A1'].load_liquid(liquid=bsa_reag_a, volume=20000)  
    reservoir['A3'].load_liquid(liquid=bsa_reag_b, volume=20000)  
    reservoir['A5'].load_liquid(liquid=bsa_reag_c, volume=20000)  
    reservoir['A7'].load_liquid(liquid=excess_lysis, volume=15000) 

    # Load pipettes
    p50_multi = protocol.load_instrument('flex_8channel_50', 'left') 
    p1000_multi = protocol.load_instrument('flex_8channel_1000', 'right') 

    #Configure the p1000 pipette to use all channels
    p1000_multi.configure_nozzle_layout(style=ALL, tip_racks=[tips_200])

    #Start recording the video
    video_process = subprocess.Popen(["python3", "/var/lib/jupyter/notebooks/record_video_western.py"])

    # assign sample locations dynamically
    sample_locations = []
    for i in range(protocol.params.num_samples):
        if i < 6:          # B1 to B6
            sample_locations.append(f'B{i + 1}')
        elif i < 12:       # C1 to C6
            sample_locations.append(f'C{i - 5}')
        elif i < 18:       # D1 to D6
            sample_locations.append(f'D{i - 11}')
        elif i < 24:              # E1 to E6 (i = 18..23)
            sample_locations.append(f'E{i - 17}')
        else:
            break  # Stop if we exceed the number of available rows/columns

    # Predefined list of letters A-H
    row = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']

    # Create a list of rows that repeats based on num (destination rows A–H cycling)
    rows = [row[i % len(row)] for i in range(protocol.params.num_samples)]

    # Create a dynamic sample map based on the assigned sample locations
    sample_map = list(map(lambda i, j: (i, j), rows, sample_locations))
    
    # ---------------- Normalizing BCA Assay ----------------
    protocol.comment("Place BCA assay absorbance data in /var/lib/jupyter/notebooks/TWH, load new deep well plate into flex B2 (where BCA plate was), and new tube rack into A2 (with excess lysis buffer in A1 and empty falcon in A2)")

    # Set P1000 to use single nozzle and mix reducing agent with loading buffer
    p1000_multi.configure_nozzle_layout(style=SINGLE, start="A1", tip_racks=[tips_200]) 

    # Add loading buffer to empty eppendorf
    p1000_multi.distribute((protocol.params.final_volume/3) * (protocol.params.num_samples + 3),
                    temp_adapter['A2'],
                    temp_adapter['A4'],
                    rate=speed,
                    mix_before=(1, 10),
                    disposal_vol=5,
                    new_tip='always')

    #Configure the p1000 and p50 pipettes to use single tip NOTE: this resets the pipettes tip racks!
    p50_multi.configure_nozzle_layout(style=SINGLE, start="A1", tip_racks=[partial_50]) #, 

    # Add Reducing agent to loading buffer
    p50_multi.distribute(((protocol.params.final_volume/3) * (protocol.params.num_samples + 3))/9,
                    temp_adapter['A3'],
                    temp_adapter['A4'],
                    rate=speed,
                    mix_before=(1, 10),
                    disposal_vol=1,
                    new_tip='always')

     # Define the directory path
    directory = Path("/var/lib/jupyter/notebooks/Data/")

    # Get today's date in YYMMDD format
    today_date = datetime.date.today().strftime("%y%m%d")

    # For debugging, change the file from wait_for_file.py to wait_for_file_debug.py
    find_file = subprocess.Popen(['python3',"/var/lib/jupyter/notebooks/wait_for_file.py"],stdout=subprocess.PIPE,
        text=True)
    stdout, stderr = find_file.communicate()

    if stderr:
        raise ValueError(f"Error while waiting for file: {stderr}")

    # Extract the file path from the output
    file_path = stdout.splitlines()[1]
    if not file_path:
        raise ValueError("No file path returned by wait_for_file.py")

    protocol.comment(f"Successfully loaded: {file_path}")
    # Read the data file
    df = pd.read_excel(file_path, header=5, nrows=8, usecols="C:N")

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

    unknown_samples = final_df.iloc[8:8 + protocol.params.num_samples]
    unknown_samples['Mean Absorbance'] = unknown_samples[['Replicate 1', 'Replicate 2', 'Replicate 3']].mean(axis=1)
    unknown_samples['Protein Concentration (mg/mL)'] = (unknown_samples['Mean Absorbance'] - intercept) / slope
    unknown_samples['Sample Volume (µL)'] = (target_concentration * protocol.params.final_volume) / unknown_samples['Protein Concentration (mg/mL)']
    unknown_samples['Diluent Volume (µL)'] = protocol.params.final_volume - unknown_samples['Sample Volume (µL)']

    # Volume check
    if any(unknown_samples['Sample Volume (µL)'] > protocol.params.final_volume):
        protocol.comment("One or more samples exceed the maximum allowed volume for dilution.")
        raise Exception("Aborting protocol: at least one sample volume exceeds the final volume threshold.")
    unknown_samples.loc[unknown_samples['Sample Volume (µL)'] > protocol.params.final_volume, ['Sample Volume (µL)', 'Diluent Volume (µL)']] = [protocol.params.final_volume, 0]
    protocol.comment("\nNormalized Unknown Samples (to 1 mg/mL in 500 µL):")
    normalized_samples = unknown_samples[['Sample', 'Protein Concentration (mg/mL)', 'Sample Volume (µL)','Diluent Volume (µL)']].reset_index().drop(columns='index')

    # Write the output and image of data plot to the instrument jupyter notebook directory
    filename = f"Protocol_output_{today_date}.csv"
    output_file_destination_path = directory.joinpath(filename)
    normalized_samples.to_csv(output_file_destination_path)
    print(unknown_samples[['Sample', 'Protein Concentration (mg/mL)', 'Sample Volume (µL)','Diluent Volume (µL)']])

    # Dilute sample in lysis buffer to 1 mg/ml on deep well plate
    rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    destination_wells  = [f'{rows[i % 8]}{(i // 8)+ 1}' for i in range(len(normalized_samples))]

    for i, row in normalized_samples.iterrows():
        source_well = sample_locations[i]
        normalized_volume = row['Sample Volume (µL)']
        diluent_volume = protocol.params.final_volume - normalized_volume
        destination_well = destination_wells[i]
        p50_multi.transfer(normalized_volume, 
                    temp_adapter[source_well], 
                    plate3[destination_well], 
                    rate=0.5, 
                    new_tip='once')
        p50_multi.transfer(diluent_volume, 
                    reservoir['A7'], 
                    plate3[destination_well], 
                    rate=0.5, 
                    new_tip='once')
    

    # Add loading buffer
    p50_multi.distribute(protocol.params.final_volume/3,
                    temp_adapter['A4'],
                    [plate3[well] for well in destination_wells],
                    rate=speed,
                    mix_after=(3, 10),
                    disposal_vol=1,
                    new_tip='always')

    # Step 3: Run thermocycling conditions
    thermocycler.close_lid()
    thermocycler.set_lid_temperature(70)
    protocol.comment('Running thermocycler for 10 minutes')
    thermocycler.set_block_temperature(70,block_max_volume=50, hold_time_minutes=10)
    thermocycler.set_block_temperature(10)  # Hold at 4°C
    # Stop video recording after the main task is completed
    video_process.terminate()