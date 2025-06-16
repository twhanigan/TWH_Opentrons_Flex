from opentrons import protocol_api
from opentrons.protocol_api import SINGLE, ALL
import pandas as pd
import numpy as np
import subprocess
from pathlib import Path
import datetime
import time

metadata = {
    'protocolName': 'Gel-based Chemical Proteomics 06132025',
    'author': 'Assistant',
    'description': 'Serial dilution of BSA standard and sample processing. This includes cooling samples to 4c, heating plate to 37c with shaking and recording a video of the whole process. Place BSA Standard in A1, Lysis buffer in A2, change the number of samples and place samples in row B starting at B1. MINIMUM Sample volumen in eppendorf tubes is 40 uL. '
}

requirements = {
    "robotType": "Flex",
    "apiLevel": "2.21"
}

def run(protocol: protocol_api.ProtocolContext):
    protocol.comment(
        "Place BSA Standard in A1, Lysis buffer in A2, tbta in A3, biotin in A4, cuso4 in A5, tcep in A6 and samples in row B")
    protocol.comment("Running the BCA assay")

    target_concentration = 1.5 #mg/ml
    final_volume = 50 # uL
    num_samples = 10 # change this to the number of samples you need to run. The maximum is 18.
    num_rows = 8  # A-H
    num_replicates = 3  # the number of replicates
    speed= 0.3 #Speed of pipetting NP40 lysis buffer=0.35, 2M Urea in EPPS=0.3

    #Start recording the video
    video_process = subprocess.Popen(["python3", "/var/lib/jupyter/notebooks/record_video.py"])

    # Load modules
    heater_shaker = protocol.load_module('heaterShakerModuleV1', 'D1')
    thermocycler = protocol.load_module('thermocyclerModuleV2')
    temp_module = protocol.load_module('temperature module gen2', 'C1')
    mag_block = protocol.load_module('magneticBlockV1', 'D2')
    chute = protocol.load_waste_chute()

    # Load adapters
    hs_adapter = heater_shaker.load_adapter('opentrons_universal_flat_adapter')
    temp_adapter = temp_module.load_labware('opentrons_24_aluminumblock_nest_1.5ml_screwcap')

    #set the heater_shaker temp to 60C
    heater_shaker.open_labware_latch()

    #set the temp module to 0c
    temp_module.set_temperature(celsius=10)
    
    # Load labware
    partial_50 = protocol.load_labware(load_name="opentrons_flex_96_filtertiprack_50ul",location="B3")
    tips_200 = protocol.load_labware(load_name="opentrons_flex_96_filtertiprack_200ul",location="A3")
    tips_1000 = protocol.load_labware('opentrons_flex_96_filtertiprack_1000ul', 'C4')
    plate1 = protocol.load_labware('opentrons_96_wellplate_200ul_pcr_full_skirt', 'A2') 
    plate2 = protocol.load_labware('corning_96_wellplate_360ul_flat', location='B2') #on heatshaker
    plate3 = protocol.load_labware('opentrons_96_wellplate_200ul_pcr_full_skirt', location='A4')  # New deep well plate for final samples
    reservoir = protocol.load_labware('nest_12_reservoir_15ml', 'C2')
    
    # Liquid definitions
    bsa_standard = protocol.define_liquid(name = 'BSA Standard', display_color="#704848",)
    bsa_reag_a = protocol.define_liquid(name = 'Reagent A', display_color="#704300",)
    bsa_reag_b = protocol.define_liquid(name = 'Reagent B', display_color="#704900",)
    bsa_reag_c = protocol.define_liquid(name = 'Reagent C', display_color="#701100",)
    excess_lysis = protocol.define_liquid(name='Excess Lysis Buffer', display_color="#FFC0CB")  # Pink
    copper_sulfate = protocol.define_liquid(name = 'copper_sulfate', display_color="#704848",)
    rhodamine_azide = protocol.define_liquid(name = 'rhodamine_azide', display_color="#704300",)
    tcep_click = protocol.define_liquid(name = 'tcep_click', display_color="#704900",)
    tbta = protocol.define_liquid(name = 'tbta', display_color="#701100",)
    loading_buffer = protocol.define_liquid(name = 'Loading Buffer', display_color="#4169E1",)
    sample_liquids = [protocol.define_liquid(name = f'Sample {i + 1}', display_color="#FFA000",) for i in range(num_samples)]

    # Reservoir assignments for washes and digestion
    reservoir['A1'].load_liquid(liquid=bsa_reag_a, volume=20000)  
    reservoir['A3'].load_liquid(liquid=bsa_reag_b, volume=20000)  
    reservoir['A5'].load_liquid(liquid=bsa_reag_c, volume=20000)  
    reservoir['A7'].load_liquid(liquid=excess_lysis, volume=15000) 
    reservoir['A9'].load_liquid(liquid=loading_buffer, volume=1000)  

    temp_adapter['A1'].load_liquid(liquid=bsa_standard, volume = 1500)
    temp_adapter['D1'].load_liquid(liquid=copper_sulfate, volume=1500) #click
    temp_adapter['D2'].load_liquid(liquid=rhodamine_azide, volume=1500) #click
    temp_adapter['D3'].load_liquid(liquid=tcep_click, volume=1500) #click
    temp_adapter['D4'].load_liquid(liquid=tbta, volume=1500) #click

    # Load pipettes
    p50_multi = protocol.load_instrument('flex_8channel_50', 'left') 
    p1000_multi = protocol.load_instrument('flex_8channel_1000', 'right') 

    #Configure the p1000 pipette to use all channels
    p1000_multi.configure_nozzle_layout(style=ALL, tip_racks=[tips_200])

    # Steps 1: Add lysis buffer to column 1 of plate1. 
    p1000_multi.distribute(50, 
         reservoir['A7'],
         plate1['A1'],
         rate = speed,
         mix_before=(1, 50),
         delay = 2,
         new_tip='once')

    #Step 3: Configure the p50 pipette to use single tip NOTE: this resets the pipettes tip racks! it doesn't
    p50_multi.configure_nozzle_layout(style=SINGLE, start="A1",tip_racks=[partial_50])

    # Step 4: Transfer BSA standard (20 mg/ml) to first well of column 1
    p50_multi.transfer(50,
        temp_adapter['A1'],
        plate1['A1'],
        rate = 0.35,
        delay = 2,
        mix_after=(3, 40),
        new_tip='once')

    # Step 5: Perform serial dilution down column 1
    rows = ['A','B', 'C', 'D', 'E', 'F', 'G']
    p50_multi.pick_up_tip()
    for source, dest in zip(rows[:-1], rows[1:]):
        p50_multi.transfer(50,
                         plate1[f'{source}1'],
                         plate1[f'{dest}1'],
                         rate = 0.5,
                         mix_after=(3, 40),
                         new_tip='never', 
                         disposal_vol=0)

    # Step 6: remove excess standard from well G
    p50_multi.aspirate(50,plate1['G1'])
    p50_multi.drop_tip()

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

    # Predefined list of letters A-H
    row = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']

    # Create a list of rows that repeats based on num_samples
    rows = [row[i % len(row)] for i in range(num_samples)]

    # Create a dynamic sample map based on the assigned sample locations
    sample_map = list(map(lambda i,j :(i,j), rows, sample_locations))

# Iterate over the sample_map list
    for index, (row, tube) in enumerate(sample_map):
        if index < 8:
            base_column = 4 + (index // 8)  # This will determine the starting column for each row
        elif index< 16:
            base_column = 6 + (index // 8)
        else:
            base_column = 8 + (index //8)

        # Prepare destination wells
        destination_wells = [f'{row}{base_column + (i % 3)}' for i in range(3)]  # Generate wells like A4, A5, A6 or B4, B5, B6, etc.
        
        #Transfer the samples onto plate 2
        p50_multi.distribute(5,
                        temp_adapter[tube],
                        [plate2[i].bottom(z=0.3) for i in destination_wells],
                        rate = speed,
                        mix_before=(1, 10),
                        disposal_vol=5)  # Distributing to three consecutive columns

    #Step 9: Load the p50 with full tip rack (don't need to)
    p50_multi.configure_nozzle_layout(style=ALL, tip_racks=[partial_50]) #, 

    #Step 10: Pipette triplicate of controls from plate1 column 1 to plate2 columns 1,2,3 
    p50_multi.distribute(5, 
                        plate1['A1'], 
                        [plate2[f'A{i}'].bottom(z=0.1) for i in range(1, 4)],
                        rate= speed,
                        mix_before=(1, 10),
                        disposal_vol=5)

    
    protocol.move_labware(labware=tips_1000, new_location='C3', use_gripper=True)
    #Step 12: Load the p1000 with full tip rack (don't need to)
    p1000_multi.configure_nozzle_layout(style=ALL, tip_racks=[tips_1000]) #,

    # Step 13: Add reagent A
    p1000_multi.distribute(50,
                        reservoir['A1'],
                        plate2.wells(),
                        new_tip='once',
                        disposal_vol=50)

    # Step 14: Add reagent B
    p1000_multi.distribute(48,
                        reservoir['A3'],
                        plate2.wells(),
                        new_tip='once',
                        disposal_vol=50)

    # Step 15: Add reagent c
    p50_multi.distribute(2,
                        reservoir['A5'],
                        plate2.wells(),
                        new_tip='once',
                        rate = speed,
                        mix_after=(2, 10),
                        disposal_vol=5)

    #Step 16: move plate 2 to the heater shaker and incubate at 37c
    protocol.move_labware(labware=plate2, new_location=hs_adapter, use_gripper=True)
    heater_shaker.set_and_wait_for_temperature(50)
    heater_shaker.close_labware_latch()
    heater_shaker.set_and_wait_for_shake_speed(500)
    protocol.delay(minutes=5)

    #Step 17 deactivate heater shaker and temp modules
    heater_shaker.deactivate_shaker()
    heater_shaker.deactivate_heater()
    heater_shaker.open_labware_latch()

    # ---------------- Normalizing BCA Assay ----------------
    # Tell the user to load BCA assay data
    protocol.comment("Place BCA assay absorbance data in /var/lib/jupyter/notebooks/Data, load new deep well plate into flex B2 (where BCA plate was), and new tube rack into A2 (with excess lysis buffer in A1 and empty falcon in A2)")

    # Pause the protocol until the user loads the file to /var/lib/jupyter/notebooks
    protocol.pause()

    # Tell the robot that new labware will be placed onto the deck
    protocol.move_labware(labware=plate1, new_location=protocol_api.OFF_DECK)
    protocol.move_labware(labware=plate2, new_location=protocol_api.OFF_DECK)
    protocol.move_labware(labware=plate3, new_location="B2", use_gripper=True)
    protocol.move_labware(labware=partial_50, new_location='B4', use_gripper=True)

    #Configure the p1000 pipette to use single tip NOTE: this resets the pipettes tip racks!
    p1000_multi.configure_nozzle_layout(style=SINGLE, start="A1",tip_racks=[tips_200])

    # Define the directory path
    directory = Path("/var/lib/jupyter/notebooks/Data/")

    # Get today's date in YYMMDD format
    today_date = datetime.date.today().strftime("%y%m%d")

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

    unknown_samples = final_df.iloc[8:8 + num_samples]
    unknown_samples['Mean Absorbance'] = unknown_samples[['Replicate 1', 'Replicate 2', 'Replicate 3']].mean(axis=1)
    unknown_samples['Protein Concentration (mg/mL)'] = (unknown_samples['Mean Absorbance'] - intercept) / slope


    unknown_samples['Sample Volume (mL)'] = (target_concentration * final_volume) / unknown_samples['Protein Concentration (mg/mL)']
    unknown_samples['Diluent Volume (mL)'] = final_volume - unknown_samples['Sample Volume (mL)']
    unknown_samples.loc[unknown_samples['Sample Volume (mL)'] > final_volume, ['Sample Volume (mL)', 'Diluent Volume (mL)']] = [final_volume, 0]
    protocol.comment("\nNormalized Unknown Samples (to 1 mg/mL in {final_volume} µL):")
    print(unknown_samples[['Sample', 'Protein Concentration (mg/mL)', 'Sample Volume (mL)', 'Diluent Volume (mL)']])

    normalized_samples = unknown_samples[['Sample', 'Protein Concentration (mg/mL)', 'Sample Volume (mL)', 'Diluent Volume (mL)']].reset_index().drop(columns='index')
    rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    destination_wells  = [f'{rows[i % 8]}{(i // 8)+ 1}' for i in range(len(normalized_samples))]

    for i, row in normalized_samples.iterrows():
        source_well = sample_locations[i]
        normalized_volume = row['Sample Volume (mL)']
        diluent_volume = final_volume - normalized_volume
        destination_well = destination_wells[i]
        p1000_multi.transfer(normalized_volume, temp_adapter[source_well], plate3[destination_well], rate=0.5, new_tip='once')
        p1000_multi.transfer(diluent_volume, reservoir['A7'], plate3[destination_well], rate=0.5, new_tip='once')

    # ---------------- Click Reaction ----------------
    protocol.comment("Running click reaction")
    protocol.move_labware(labware=partial_50, new_location='B3', use_gripper=True)
    protocol.move_labware(labware=tips_1000, new_location='C4', use_gripper=True)
    p50_multi.configure_nozzle_layout(style=SINGLE, start="A1", tip_racks=[partial_50]) #,
    
    #Pipette rhodamine azide (d2), tbta (d4), cuso4 (d1), and tcep (d3)
    p50_multi.transfer(1*(num_samples*1.5), 
                            temp_adapter['D2'], 
                            temp_adapter['D6'],
                            rate=speed,
                            mix_before=(1,10), 
                            new_tip='once')

    p50_multi.transfer(6*(num_samples*1.5), 
                            temp_adapter['D4'], 
                            temp_adapter['D6'],
                            mix_before=(1,10),
                            rate=speed,
                            delay=3, 
                            new_tip='once')

    p50_multi.transfer(2*(num_samples*1.5), 
                            temp_adapter['D1'], 
                            temp_adapter['D6'], 
                            new_tip='once')

    p50_multi.transfer(2*(num_samples*1.5), 
                            temp_adapter['D3'], 
                            temp_adapter['D6'], 
                            new_tip='once')
    
    # Pipette the click reaction premix
    p50_multi.distribute(6, 
                            temp_adapter['D6'], 
                            [plate3[i] for i in destination_wells],
                            rate=speed-0.05,
                            mix_before=(3, 30),
                            mix_after=(3, 30), 
                            new_tip='always')

    # Step 11: shake the sample plate for click reaction
    protocol.move_labware(labware=plate3, new_location=hs_adapter, use_gripper=True)
    heater_shaker.close_labware_latch()
    heater_shaker.set_and_wait_for_shake_speed(1000)
    protocol.delay(minutes=60)
    heater_shaker.deactivate_shaker()

    # Add the loading buffer and move to the thermocylcer to seal and store.
    p50_multi.configure_nozzle_layout(style=ALL, tip_racks=[partial_50])
    p50_multi.distribute(34, 
                            reservoir['A9'], 
                            plate3[destination_wells], 
                            mix_after=(3, 40), 
                            new_tip='always')

    heater_shaker.open_labware_latch()
    thermocycler.open_lid()
    protocol.move_labware(labware=plate3, new_location=thermocycler, use_gripper=True)
    thermocycler.close_lid()
    thermocycler.set_block_temperature(4)  # Hold at 4°C
