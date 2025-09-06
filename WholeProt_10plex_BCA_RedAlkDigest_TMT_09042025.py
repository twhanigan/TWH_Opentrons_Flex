from opentrons import protocol_api
from opentrons.protocol_api import SINGLE, ALL
import pandas as pd
import numpy as np
#import matplotlib.pyplot as plt
import subprocess
from pathlib import Path
import datetime
import time

metadata = {
    'protocolName': 'TMT-labeled whole proteome sample preparation 09042025',
    'author': 'Assistant',
    'description': 'I changed the tips in this one. BCA assay and sample protein normalization (<20 Samples), reduction/alkylation and SP3 sample cleanup, trypsin digestion and on bead TMT reaction.'
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
        default=10,
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

    parameters.add_str(
        variable_name="tmt_row",
        display_name="TMT row",
        description="Row on the TMT plate to use (A–H)",
        default="A",
        choices=[{"display_name": r, "value": r} for r in "ABCDEFGH"]
    )

    parameters.add_float(
        variable_name="target_concentration",
        display_name="Target protein concentration",
        description="Concentration of normalized lysate in click reaction",
        default=1,
        minimum=0.5,
        maximum=2.5,
        unit="µg/µL"
    )

    parameters.add_int(
        variable_name="final_volume",
        display_name="final volume lysate",
        description="Volume to normalize µg of protein in",
        default=50,
        minimum=20,
        maximum=1000,
        unit="µL"
    )

def run(protocol: protocol_api.ProtocolContext):
    # ---------------- BCA Assay for Preprepared Lysates ----------------
    protocol.comment("Running BCA. KEEP ON ROCKING!")
    
    #Speed of pipetting NP40 lysis buffer=0.35, 2M Urea in EPPS=0.3
    speed= 0.3
    num_rows = 8  # A-H

    #Start recording the video
    video_process = subprocess.Popen(["python3", "/var/lib/jupyter/notebooks/record_video_chemprot.py"])

    # Load modules
    heater_shaker = protocol.load_module('heaterShakerModuleV1', 'D1')
    thermocycler = protocol.load_module('thermocyclerModuleV2')
    temp_module = protocol.load_module('temperature module gen2', 'C1')
    mag_block = protocol.load_module('magneticBlockV1', 'D2')
    chute = protocol.load_waste_chute()

    # Load adapters
    temp_adapter = temp_module.load_labware('opentrons_24_aluminumblock_nest_1.5ml_screwcap')
    epp_rack = protocol.load_labware('opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap', location=protocol_api.OFF_DECK)

    #set the heater_shaker temp to 60C
    thermocycler.open_lid()
    heater_shaker.open_labware_latch()

    #set the temp module to 0c
    temp_module.set_temperature(celsius=10)
    
    # Load labware
    tips_50 = protocol.load_labware(load_name="opentrons_flex_96_filtertiprack_50ul",location="A3")
    tips_200 = protocol.load_labware(load_name='opentrons_flex_96_filtertiprack_200ul', location='B3')
    tips_1000 = protocol.load_labware(load_name="opentrons_flex_96_filtertiprack_1000ul",location="C4")
    tips_1000_2 = protocol.load_labware(load_name="opentrons_flex_96_filtertiprack_1000ul",location="D4")
    plate1 = protocol.load_labware('opentrons_96_wellplate_200ul_pcr_full_skirt', 'A2')
    plate2 = protocol.load_labware('corning_96_wellplate_360ul_flat', 'B2')
    plate3 = protocol.load_labware('nest_96_wellplate_2ml_deep', location=protocol_api.OFF_DECK) 
    plate4 = protocol.load_labware('opentrons_96_wellplate_200ul_pcr_full_skirt', location=protocol_api.OFF_DECK)
    tmt_plate = protocol.load_labware('opentrons_96_wellplate_200ul_pcr_full_skirt', location=protocol_api.OFF_DECK)
    reservoir = protocol.load_labware('nest_12_reservoir_15ml', 'C2')
    
    # load the liquids to the tube racks and reservoirs
    bsa_standard = protocol.define_liquid(name='BSA Standard', display_color="#FF00FF")  # Brown
    bsa_reag_a = protocol.define_liquid(name = 'Reagent A', display_color="#00FFFF")
    bsa_reag_b = protocol.define_liquid(name = 'Reagent B', display_color="#FFFF00")
    bsa_reag_c = protocol.define_liquid(name = 'Reagent C', display_color="#FF3300")
    SP3_Beads = protocol.define_liquid(name='SP3 Beads', display_color="#00FF00")  # Red
    K2CO3 = protocol.define_liquid(name='K2CO3', display_color="#FF6600")  # Purple
    tcep = protocol.define_liquid(name='TCEP', display_color="#0066FF")  # Orange
    empty_2mL = protocol.define_liquid(name='empty_2mL', display_color="#1685A4")  # SteelBlue
    IAA =  protocol.define_liquid(name='IAA', display_color="#AA00FF")  # ?
    excess_lysis = protocol.define_liquid(name='Excess Lysis Buffer', display_color="#FF0099")  # Pink
    epps_urea = protocol.define_liquid(name='2M Urea in EPPS', display_color="#00FF99")  # SaddleBrown
    epps = protocol.define_liquid(name='2M Urea in EPPS', display_color="#00FF99")  # SaddleBrown
    abs_ethanol = protocol.define_liquid(name='100% Ethanol', display_color="#4682B4")  # SteelBlue
    ethanol = protocol.define_liquid(name='80% Ethanol', display_color="#4682B4")  # SteelBlue
    ACN = protocol.define_liquid(name='Acetonitrile', display_color="#4682C3")  # SteelBlue
    trypsin = protocol.define_liquid(name='Trypsin in EPPS', display_color="#9900FF")  # Gold
    cacl2 = protocol.define_liquid(name='CaCl2', display_color="#FF3300")  # LimeGreen
    hydroxylamine = protocol.define_liquid(name='hydroxylamine', display_color="#8A2BE2")   # Blue Violet
    sample_liquids = [protocol.define_liquid(name = f'Sample {i + 1}', display_color="#FFA000",) for i in range(protocol.params.num_samples)]

    # Reservoir assignments for washes and digestion
    reservoir['A1'].load_liquid(liquid=bsa_reag_a, volume=15000)  
    reservoir['A2'].load_liquid(liquid=bsa_reag_b, volume=15000)  
    reservoir['A3'].load_liquid(liquid=bsa_reag_c, volume=15000)  
    reservoir['A4'].load_liquid(liquid=excess_lysis, volume=15000) 
    reservoir['A6'].load_liquid(liquid=ACN, volume=15000)
    reservoir['A8'].load_liquid(liquid=abs_ethanol,volume=15000)
    reservoir['A9'].load_liquid(liquid=ethanol, volume=15000)  # 80% Ethanol
    reservoir['A10'].load_liquid(liquid=ethanol, volume=15000)  # 80% Ethanol for washes
    reservoir['A11'].load_liquid(liquid=epps, volume = 15000)   # 200 mM EPPS
    reservoir['A12'].load_liquid(liquid=epps_urea, volume=15000)  # 2M Urea in EPPS

    # Trypsin and CaCl2 for digestion
    temp_adapter['A1'].load_liquid(liquid=bsa_standard, volume=1000)  # 20 mg/ml BSA standard
    temp_adapter['A2'].load_liquid(liquid=SP3_Beads, volume=1000)  # SP3 beads for desalt
    temp_adapter['A3'].load_liquid(liquid=tcep, volume=1000)  # reduce
    temp_adapter['A4'].load_liquid(liquid=K2CO3, volume=1000)  # reduce
    temp_adapter['A5'].load_liquid(liquid=IAA, volume=1000)  # alkylate
    temp_adapter['A6'].load_liquid(liquid=empty_2mL, volume=0)  # empty tube to combine TCEP/K2CO3
    temp_adapter['D4'].load_liquid(liquid=hydroxylamine, volume=200) # Quench TMT reaction
    temp_adapter['D5'].load_liquid(liquid=cacl2, volume=500)  # CaCl2 for digestion
    temp_adapter['D6'].load_liquid(liquid=trypsin, volume=2000)  # Trypsin in EPPS for digestion

    # Load pipettes
    p50_multi = protocol.load_instrument('flex_8channel_50', 'left') 
    p1000_multi = protocol.load_instrument('flex_8channel_1000', 'right')

    #Configure the p1000 pipette to use single tip NOTE: this resets the pipettes tip racks!
    p1000_multi.configure_nozzle_layout(style=ALL, tip_racks=[tips_200])

    # Steps 1: Add lysis buffer to column 1 of plate1. 
    p1000_multi.transfer(50, 
         reservoir['A4'],
         plate1[f'A{protocol.params.standards_col}'],
         rate = speed,
         delay = 2,
         new_tip='once',
         blow_out=True)

    # Step 2: move the 200uL tips to D4 and then the 50 uL partial tips to B3
    protocol.move_labware(labware=tips_200, new_location="B4", use_gripper=True)

    #Step 3: Configure the p50 pipette to use single tip NOTE: this resets the pipettes tip racks!
    p50_multi.configure_nozzle_layout(style=SINGLE, start="A1",tip_racks=[tips_50])

    # Step 4: Transfer BSA standard (20 mg/ml) to first well of column 1
    p50_multi.transfer(50,
        temp_adapter['A1'],
        plate1[f'A{protocol.params.standards_col}'],
        rate = 0.35,
        delay = 2,
        mix_after=(3, 40),
        new_tip='once')

    # Step 5: Perform serial dilution down column 1
    rows = ['A','B', 'C', 'D', 'E', 'F', 'G']
    p50_multi.pick_up_tip()
    for source, dest in zip(rows[:-1], rows[1:]):
        p50_multi.transfer(50,
                         plate1[f'{source}{protocol.params.standards_col}'],
                         plate1[f'{dest}{protocol.params.standards_col}'],
                         rate = 0.5,
                         mix_after=(3, 40),
                         new_tip='never', 
                         disposal_vol=0)

    # Step 6: remove excess standard from well G
    p50_multi.aspirate(50,plate1[f'G{protocol.params.standards_col}'])
    p50_multi.drop_tip()

    # assign sample locations dynamically
    sample_locations = []
    for i in range(protocol.params.num_samples):
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
    rows = [row[i % len(row)] for i in range(protocol.params.num_samples)]


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
                        [plate2[i].bottom(z=0.1) for i in destination_wells],
                        rate = speed,
                        mix_before=(1, 10),
                        disposal_vol=5)  # Distributing to three consecutive columns

    #Step 9: Load the p50 with full tip rack
    p50_multi.configure_nozzle_layout(style=ALL, tip_racks=[tips_50]) #, 

    #Step 10: Pipette triplicate of controls from plate1 column 1 to plate2 columns 1,2,3 
    p50_multi.distribute(5, 
                        plate1[f'A{protocol.params.standards_col}'], 
                        [plate2[f'A{i}'].bottom(z=0.1) for i in range(1, 4)],
                        rate= speed,
                        mix_before=(1, 10),
                        disposal_vol=5)

    # Step 11: move the 50 uL partial tips to A4 and the 200uL complete tips to B3
    protocol.move_labware(labware=tips_1000, new_location="B3", use_gripper=True)

    #Step 12: Load the p1000 with full tip rack
    p1000_multi.configure_nozzle_layout(style=ALL, tip_racks=[tips_1000]) #,

    # Step 13: Add reagent A
    p1000_multi.distribute(50,
                        reservoir['A1'],
                        plate2.wells(),
                        new_tip='once',
                        disposal_vol=50)

    # Step 14: Add reagent B
    p1000_multi.distribute(48,
                        reservoir['A2'],
                        plate2.wells(),
                        new_tip='once',
                        disposal_vol=50)

    # Step 15: Add reagent c
    p50_multi.distribute(2,
                        reservoir['A3'],
                        plate2.wells(),
                        new_tip='once',
                        rate = speed,
                        mix_after=(2, 10),
                        disposal_vol=5)

    #Step 16: move plate 2 to the heater shaker and incubate at 37c
    heater_shaker.open_labware_latch()
    protocol.move_labware(labware=plate2, new_location=heater_shaker,use_gripper=True)
    heater_shaker.set_and_wait_for_temperature(40)
    heater_shaker.close_labware_latch()
    heater_shaker.set_and_wait_for_shake_speed(500)
    protocol.delay(minutes=10)

    #Step 17 deactivate heater shaker and temp modules
    heater_shaker.deactivate_shaker()
    heater_shaker.deactivate_heater()
    heater_shaker.open_labware_latch()
    protocol.move_labware(labware=plate2, new_location='B2',use_gripper=True)


    # ---------------- Normalizing BCA Assay ----------------
    protocol.comment("Place BCA assay absorbance data in /var/lib/jupyter/notebooks/TWH, load new deep well plate into hs_adapter, and new tube rack into A2 (with click mix, beads and TCEP/IAA)")
    
    # Pause the protocol until the user loads the file to /var/lib/jupyter/notebooks
    protocol.pause()

    # Tell the robot that new labware will be placed onto the deck
    protocol.move_labware(labware=plate1, new_location=protocol_api.OFF_DECK)
    protocol.move_labware(labware=plate2, new_location=protocol_api.OFF_DECK)
    protocol.move_labware(labware=plate3,new_location='B2')
    protocol.move_labware(labware=plate4, new_location=thermocycler)
    protocol.move_labware(labware=tmt_plate, new_location='A2')

    #Define the directory path
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
    unknown_samples['Sample Volume (µL)'] = (protocol.params.target_concentration * protocol.params.final_volume) / unknown_samples['Protein Concentration (mg/mL)']
    unknown_samples['Diluent Volume (µL)'] = protocol.params.final_volume - unknown_samples['Sample Volume (µL)']
    unknown_samples.loc[unknown_samples['Sample Volume (µL)'] > protocol.params.final_volume, ['Sample Volume (µL)', 'Diluent Volume (µL)']] = [protocol.params.final_volume, 0]
    protocol.comment(f"\nNormalized Unknown Samples (to {protocol.params.target_concentration} mg/mL in {protocol.params.final_volume} µL):")
    summary = unknown_samples[['Sample', 'Protein Concentration (mg/mL)', 'Sample Volume (µL)', 'Diluent Volume (µL)']].to_string(index=False)
    protocol.comment(f"\nNormalized sample volumes:\n{summary}")

    normalized_samples = unknown_samples[['Sample', 'Protein Concentration (mg/mL)', 'Sample Volume (µL)', 'Diluent Volume (µL)']].reset_index().drop(columns='index')
    # Write the output and image of data plot to the instrument jupyter notebook directory
    filename = f"Protocol_output_{today_date}.csv"
    output_file_destination_path = directory.joinpath(filename)
    normalized_samples.to_csv(output_file_destination_path)
    rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    destination_wells  = [f'{rows[i % 8]}{(i // 8)+ 1}' for i in range(len(normalized_samples))]

    # Dilute sample in lysis buffer to 1 mg/ml on deep well plate
    rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    destination_wells  = [f'{rows[i % 8]}{(i // 8)+ 1}' for i in range(len(normalized_samples))]

    tiprack_locations = {
        "tips_50": "A3",  # assume starts here
        "tips_200":"B4",    # assume starts her
        "tips_1000": "B3",   # assume starts here
        "tips_1000_2": "D4"
    }
    def normalize_samples():
        pipette = None
        if protocol.params.final_volume < 50:
            if tiprack_locations["tips_1000"] != "C4":
                protocol.move_labware(tips_1000, new_location="C4", use_gripper=True)
                tiprack_locations["tips_1000"] = "C4"
            p50_multi.configure_nozzle_layout(style=SINGLE, start="A1", tip_racks=[tips_50])
            pipette = p50_multi
        elif 50 < protocol.params.final_volume < 200:
            if tiprack_locations["tips_1000"] != "C4":
                protocol.move_labware(tips_1000, new_location="C4", use_gripper=True)
                tiprack_locations["tips_1000"] = "C4"
            if tiprack_locations["tips_200"] != "B3":
                protocol.move_labware(tips_200, new_location="B3", use_gripper=True)
                tiprack_locations["tips_200"] = "B3"
            p1000_multi.configure_nozzle_layout(style=SINGLE, start="A1", tip_racks=[tips_200])
            pipette = p1000_multi
        elif 200 < protocol.params.final_volume:
            if tiprack_locations["tips_1000"] != "B3":
                protocol.move_labware(tips_1000, new_location="B3", use_gripper=True)
                tiprack_locations["tips_1000"] = "B3"
            p1000_multi.configure_nozzle_layout(style=SINGLE, start="A1", tip_racks=[tips_1000])
            pipette = p1000_multi

        for i, row in normalized_samples.iterrows():
            source_well = sample_locations[i]
            normalized_volume = row['Sample Volume (µL)']
            diluent_volume = protocol.params.final_volume - normalized_volume
            destination_well = destination_wells[i]
            pipette.transfer(normalized_volume, 
                                temp_adapter[source_well], 
                                plate3[destination_well].bottom(z=0.1),
                                rate=speed, 
                                new_tip='once')
            pipette.transfer(diluent_volume, 
                                reservoir['A4'], 
                                plate3[destination_well].bottom(z=0.1), 
                                rate=speed,
                                new_tip='once')
    normalize_samples()

    # ---------------- Reduction, Alkylation and Digestion ----------------
    #set the heater_shaker temp to 37 c for reduce
    heater_shaker.set_and_wait_for_temperature(37)

    #Configure the p1000 pipette to use single tip NOTE: this resets the pipettes tip racks!
    if tiprack_locations["tips_200"] != "B4":
        protocol.move_labware(tips_50, new_location="B4", use_gripper=True)
        tiprack_locations["tips_50"] = "B4"
    if tiprack_locations["tips_1000"] != "C4":
        protocol.move_labware(tips_1000, new_location="C4", use_gripper=True)
        tiprack_locations["tips_1000"] = "C4"
    p1000_multi.configure_nozzle_layout(style=SINGLE, start="A1",tip_racks=[tips_1000])

    #Pipette mix K2CO3, and TCEP and transfer to samples above the 200 uL liquid miniscus.
    Redu_mix = protocol.params.num_samples*5*1.2/2 # 10 mM TCEP, 30 mM K2CO2 final conc from a 200 mM/600 mM stock
    p50_multi.transfer(Redu_mix, 
                        temp_adapter['A3'], 
                        temp_adapter['A6'], 
                        rate=speed,
                        mix_before=(1,50),
                        new_tip='once')

    p50_multi.transfer(Redu_mix, 
                        temp_adapter['A4'], 
                        temp_adapter['A6'], 
                        rate=speed,
                        mix_before=(1,50),
                        new_tip='always')

    p50_multi.distribute(5, 
                        temp_adapter['A6'], 
                        [plate3[i].bottom(z=1) for i in destination_wells],
                        mix_before=(1,50), 
                        rate=speed,
                        new_tip='once')
    
    # Reduce for 30 minutes
    protocol.move_labware(labware=plate3, new_location=heater_shaker, use_gripper=True)
    heater_shaker.close_labware_latch()
    heater_shaker.set_and_wait_for_shake_speed(1000)
    protocol.delay(minutes=30)
    heater_shaker.deactivate_shaker()
    heater_shaker.open_labware_latch()

    # add IAA and alkylate for 30 minutes
    protocol.move_labware(labware=plate3, new_location='B2', use_gripper=True)
    p50_multi.distribute(7, 
                            temp_adapter['A5'], 
                            [plate3[i].bottom(z=1) for i in destination_wells],
                            rate=speed,
                            mix_before=(1,50),
                            new_tip='once')

    protocol.move_labware(labware=plate3, new_location=heater_shaker, use_gripper=True)
    heater_shaker.close_labware_latch()
    heater_shaker.set_and_wait_for_shake_speed(1000)
    protocol.delay(minutes=30)
    heater_shaker.deactivate_shaker()

    # ---------------- Protein Cleanup by SP3 ----------------
    protocol.comment("Desalting")

    # mix the sp3 beads to homogenize. For a 10-plex sample, you need atleast 600 uL of beads.
    protocol.move_labware(tips_1000, new_location="B3", use_gripper=True)
    tiprack_locations["tips_1000"] = "B3"
    p1000_multi.configure_nozzle_layout(style=SINGLE, start="A1", tip_racks=[tips_1000])
    p1000_multi.pick_up_tip()
    p1000_multi.mix(4, 200, temp_adapter['A2'])
    p1000_multi.drop_tip()

    # transfer the sp3 beads for protein precipitation
    p1000_multi.distribute(50, 
                            temp_adapter['A2'], 
                            [plate3[i] for i in destination_wells],
                            rate=speed)
    
    # incubate beads with shaking for 5 minutes
    heater_shaker.open_labware_latch()
    protocol.move_labware(labware=plate3, new_location=heater_shaker, use_gripper=True)
    heater_shaker.close_labware_latch()
    heater_shaker.set_and_wait_for_shake_speed(1000)
    protocol.delay(minutes=5)
    heater_shaker.deactivate_shaker()
    
    #Move the full rack of p1000 tips
    protocol.move_labware(labware=tips_1000_2, new_location="C3", use_gripper=True)
    tiprack_locations["tips_1000_2"] = "C3"

    #Configure the p1000 pipette to use All tips 
    p1000_multi.configure_nozzle_layout(style=ALL, tip_racks=[tips_1000_2])

    # Determine number of full columns to fill and convert columns to top row wells
    num_full_columns = (protocol.params.num_samples + 7) // 8  # Round up to ensure all samples are covered
    destination_columns = plate3.columns()[:num_full_columns]
    destination_wells_col = [col[0] for col in destination_columns]  # Only use top row for multi-channel pipette

    # Add 100% EtOH to the sample columns and use air gap for volatility Total volume of 380 with beads/TCEP/IAA/Sample
    p1000_multi.distribute(268, 
                            reservoir['A8'], 
                            [well.top() for well in destination_wells_col], 
                            new_tip='once',
                            rate=speed,
                            air_gap=10)

    # incubate ethanol with shaking for 5 minutes
    heater_shaker.set_and_wait_for_shake_speed(1000)
    protocol.delay(minutes=5)
    heater_shaker.deactivate_shaker()
    
    # Move samples to the magnet
    heater_shaker.open_labware_latch()
    protocol.move_labware(labware=plate3, new_location=mag_block, use_gripper=True)
    protocol.delay(minutes=18) ##this is a minimal time necessary to bind the beads in that much EtOH/lysis buffer.

    # Remove the EtOH from the beads leaving 30 uL in the bottom
    p1000_multi.pick_up_tip()
    for well in destination_wells_col:
        p1000_multi.aspirate(200, well.bottom(z=2), rate=0.25)
        p1000_multi.blow_out(chute)

    # Remove remainder of EtOH and Add 80% EtOH and wash the beads three times with 200 uL 80% EtOH
    for i in range(3):
        for well in destination_wells_col:
            p1000_multi.aspirate(180, well.bottom(z=1), rate=speed-0.25)
            pipette.delay(seconds=2)
            p1000_multi.blow_out(chute)

        # add 80% EtOH, brief shake, then back to magnet
        p1000_multi.distribute(
            180,
            reservoir['A10'],
            [w.top() for w in destination_wells_col],
            rate=speed,
            new_tip='never'
        )
        protocol.move_labware(labware=plate3, new_location=heater_shaker, use_gripper=True)
        heater_shaker.close_labware_latch()

        heater_shaker.set_and_wait_for_shake_speed(1000)
        protocol.delay(seconds=60)
        heater_shaker.deactivate_shaker()
        heater_shaker.open_labware_latch()
        protocol.move_labware(labware=plate3, new_location=mag_block, use_gripper=True)
        protocol.delay(minutes=3)

        # After the final (3rd) wash, remove the last ~200 µL
        if i == 2:
            for well in destination_wells_col:
                p1000_multi.aspirate(180, well.bottom(z=0.5), rate=0.1)
                p1000_multi.blow_out(chute)
    p1000_multi.drop_tip()
    
    # Air dry
    protocol.delay(minutes=10)

    # Resuspend in 2 M urea in EPPS for digestion setup
    protocol.move_labware(labware=plate3, new_location='B2', use_gripper=True)
    p1000_multi.transfer(
        47.5,
        reservoir['A12'],
        [well for well in destination_wells_col],
        mix_before=(1, 150),
        mix_after=(3, 150),
        new_tip='always'
    )

    # Add CaCl2, trypsin in epps, and move to shaker
    p1000_multi.distribute(2.5, 
                            temp_adapter['D5'], 
                            [plate3[i].top() for i in destination_wells],
                            disposal_vol=10, 
                            new_tip='once')
    # Resuspend trypsin in 600 uL 2 M urea EPPS
    p1000_multi.transfer(50, 
                            temp_adapter['D6'], 
                            [plate3[i].top() for i in destination_wells],
                            disposal_vol=0,
                            rate=speed,
                            mix_before=(1,150),
                            mix_after=(3,500),
                            new_tip='once')

    # Digest overnight
    protocol.move_labware(labware=plate3, new_location=heater_shaker, use_gripper=True)
    heater_shaker.close_labware_latch()
    heater_shaker.set_and_wait_for_shake_speed(1000)
    protocol.delay(minutes=960)
    heater_shaker.deactivate_shaker()
    heater_shaker.open_labware_latch()
    protocol.move_labware(labware=plate3, new_location=mag_block, use_gripper=True)
    protocol.delay(minutes=3)
    protocol.comment("Samples have been digested")
   
    # ---------------- Post-Digestion TMT Reaction ----------------
    protocol.comment("Running TMT reaction")
    protocol.move_labware(labware=plate3, new_location='B2', use_gripper=True)

    # Build destination wells for plate4 starting at column 3 (C3..H3, then C4..)
    num = protocol.params.num_samples
    tmt_dilution_wells = [plate4.wells_by_name()[f"{rows[i % 8]}{(i // 8) + 3}"] for i in range(num)]

    # Use ALL nozzles for column-wise moves
    p50_multi.configure_nozzle_layout(style=ALL, tip_racks=[tips_50])

    # protein (µg/µL) * volume (µL) = µg
    total_protein = protocol.params.final_volume * protocol.params.target_concentration  # µg
    tmt_vol = 10 * (47.5 + 2.5 + 50) / total_protein   # µL needed to move 10 µg

    # Move "tmt_vol" of sample from plate3 -> plate4 and top up with EPPS to 10 µL
    for c in range(num_full_columns):
        src = plate3.columns()[c][0].bottom(0.1)  # A-row well of the column
        dst = plate4.columns()[c][0].bottom(0.1)  # A-row well of the column

        p50_multi.transfer(
            tmt_vol,
            src,
            dst,
            new_tip='always'
        )


        p50_multi.transfer(
            20 - tmt_vol,
            reservoir['A11'],
            dst,
            new_tip='always',
            mix_after=(3, 5)
        )
    # Add TMT reagents (single nozzle addressing specific wells)
    p50_multi.configure_nozzle_layout(style=SINGLE, start="A1", tip_racks=[tips_50])
    tmt_sources = [tmt_plate[f'{protocol.params.tmt_row}{i+1}'] for i in range(num)]

    for source, dest in zip(tmt_sources, tmt_dilution_wells):
        p50_multi.transfer(
            5,
            source.bottom(0.1),
            dest.bottom(0.1),
            disposal_volume=0,
            rate=speed-0.2,
            new_tip='always',
            mix_after=(3, 10)
        )

    # React
    protocol.delay(minutes=120)

    # Quench with hydroxylamine
    for dest in tmt_dilution_wells:
        p50_multi.transfer(
            1,
            temp_adapter['D4'],
            dest.bottom(1),
            rate=speed,
            mix_before=(1,5),
            mix_after=(2, 5),
            new_tip='always'   # avoid cross-contamination
        )

    thermocycler.close_lid()
    protocol.delay(minutes=5)

    video_process.terminate()

    thermocycler.set_block_temperature(4)  # hold at 4°C
