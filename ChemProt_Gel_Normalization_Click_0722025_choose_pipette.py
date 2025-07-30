from opentrons import protocol_api
from opentrons.protocol_api import SINGLE, ALL
import pandas as pd
import numpy as np
import subprocess
from pathlib import Path
import datetime
import time

metadata = {
    'protocolName': 'Gel-based Chemical Proteomics 10-sample 07232025',
    'author': 'Assistant',
    'description': 'Normalization and click reaction for Gel-based photolabeling '
}

requirements = {
    "robotType": "Flex",
    "apiLevel": "2.21"
}

def run(protocol: protocol_api.ProtocolContext):
    protocol.comment(
        "The Hanigan lab is the best lab!")
 
    target_concentration = 1.5 #mg/ml
    final_volume = 50 # uL
    num_samples = 10 # change this to the number of samples you need to run. The maximum is 18.
    num_rows = 8  # A-H
    num_replicates = 3  # the number of replicates
    speed= 0.35 #Speed of pipetting NP40 lysis buffer=0.35, 2M Urea in EPPS=0.3

    #Start recording the video
    video_process = subprocess.Popen(["python3", "/var/lib/jupyter/notebooks/record_video.py"])

    # Load modules
    heater_shaker = protocol.load_module('heaterShakerModuleV1', 'D1')
    thermocycler = protocol.load_module('thermocyclerModuleV2')
    temp_module = protocol.load_module('temperature module gen2', 'C1')
    mag_block = protocol.load_module('magneticBlockV1', 'D2')
    chute = protocol.load_waste_chute()

    # Load adapters
    #hs_adapter = heater_shaker.load_adapter('opentrons_universal_flat_adapter')
    temp_adapter = temp_module.load_labware('opentrons_24_aluminumblock_nest_1.5ml_screwcap')

    #set the heater_shaker temp to 60C
    heater_shaker.open_labware_latch()

    #set the temp module to 0c
    temp_module.set_temperature(celsius=10)
    
    # Load labware
    partial_50 = protocol.load_labware(load_name="opentrons_flex_96_filtertiprack_50ul",location="B4")
    tips_200 = protocol.load_labware(load_name="opentrons_flex_96_filtertiprack_200ul",location="A3")
    tips_1000 = protocol.load_labware('opentrons_flex_96_filtertiprack_1000ul', 'C4')
    plate3 = protocol.load_labware('opentrons_96_wellplate_200ul_pcr_full_skirt', location='B2')  # New deep well plate for final samples
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
    empty = protocol.define_liquid(name = 'empty', display_color="#FFFFFF")
    sample_liquids = [protocol.define_liquid(name = f'Sample {i + 1}', display_color="#FFA000",) for i in range(num_samples)]

    # Reservoir assignments for washes and digestion
    reservoir['A1'].load_liquid(liquid=bsa_reag_a, volume=20000)  
    reservoir['A3'].load_liquid(liquid=bsa_reag_b, volume=20000)  
    reservoir['A5'].load_liquid(liquid=bsa_reag_c, volume=20000)  
    reservoir['A7'].load_liquid(liquid=excess_lysis, volume=15000) 
    reservoir['A9'].load_liquid(liquid=loading_buffer, volume=1000)  

    temp_adapter['A1'].load_liquid(liquid=bsa_standard, volume = 1500)
    temp_adapter['A2'].load_liquid(liquid=copper_sulfate, volume=1500) #click
    temp_adapter['A3'].load_liquid(liquid=rhodamine_azide, volume=1500) #click
    temp_adapter['A4'].load_liquid(liquid=tcep_click, volume=1500) #click
    temp_adapter['A5'].load_liquid(liquid=tbta, volume=1500) #click
    temp_adapter['A6'].load_liquid(liquid=empty, volume=1500) #click

    # Load pipettes
    p50_multi = protocol.load_instrument('flex_8channel_50', 'left') 
    p1000_multi = protocol.load_instrument('flex_8channel_1000', 'right') 

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

    # ---------------- Normalizing BCA Assay ----------------
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

    def normalize_samples():
        # Create absorbance DataFrame
        initial_df = pd.DataFrame({'Well': well_names, 'Absorbance': absorbance_values})
        samples, replicate_1, replicate_2, replicate_3 = [], [], [], []
        sample_index = 1
        for col_offset in range(0, 12, 3):
            for row_offset in range(8):
                start = row_offset * 12 + col_offset
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

        # Calibration using standard samples (1–8)
        samples_1_to_8 = final_df.iloc[:8].copy()
        samples_1_to_8['Mean Absorbance'] = samples_1_to_8[['Replicate 1', 'Replicate 2', 'Replicate 3']].mean(axis=1)
        protein_concentrations = [10, 5, 2.5, 1.25, 0.625, 0.3125, 0.15625, 0]
        samples_1_to_8['Protein Concentration (mg/mL)'] = protein_concentrations
        slope, intercept = np.polyfit(
            samples_1_to_8['Protein Concentration (mg/mL)'],
            samples_1_to_8['Mean Absorbance'],
            1
        )

        # Process unknown samples
        unknown_samples = final_df.iloc[8:8 + num_samples].copy()
        unknown_samples['Mean Absorbance'] = unknown_samples[['Replicate 1', 'Replicate 2', 'Replicate 3']].mean(axis=1)
        unknown_samples['Protein Concentration (mg/mL)'] = (
            (unknown_samples['Mean Absorbance'] - intercept) / slope
        )
        unknown_samples['Sample Volume (µL)'] = (target_concentration * final_volume) / unknown_samples['Protein Concentration (mg/mL)']
        unknown_samples['Diluent Volume (µL)'] = final_volume - unknown_samples['Sample Volume (µL)']
        unknown_samples.loc[
            unknown_samples['Sample Volume (µL)'] > final_volume,
            ['Sample Volume (µL)', 'Diluent Volume (µL)']
        ] = [final_volume, 0]

        protocol.comment(f"Normalized Unknown Samples (to {target_concentration} mg/mL in {final_volume} µL):")
        summary = unknown_samples[['Sample', 'Protein Concentration (mg/mL)', 'Sample Volume (µL)', 'Diluent Volume (µL)']].to_string(index=False)
        protocol.comment(f"\nNormalized sample volumes:\n{summary}")

        # Write the output and image of data plot to the instrument jupyter notebook directory
        filename = f"Gel_chemprot_output_{today_date}.csv"
        output_file_destination_path = directory.joinpath(filename)
        normalized_samples.to_csv(output_file_destination_path)

        # Create destination well map
        rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
        destination_wells = [f'{rows[i % 8]}{(i // 8) + 1}' for i in range(len(unknown_samples))]

        # Track tip usage state
        current_tiprack = None

        def use_pipette(volume_uL):
            nonlocal current_tiprack

            if volume_uL <= 50:
                if current_tiprack != 'partial_50':
                    protocol.move_labware(partial_50, new_location='B3', use_gripper=True)
                    p50_multi.configure_nozzle_layout(style=SINGLE, start="A1", tip_racks=[partial_50])
                    current_tiprack = 'partial_50'
                return p50_multi

            elif 50 < volume_uL <= 200:
                if current_tiprack != 'tips_200':
                    protocol.move_labware(tips_200, new_location='B3', use_gripper=True)
                    p1000_multi.configure_nozzle_layout(style=SINGLE, start="A1", tip_racks=[tips_200])
                    current_tiprack = 'tips_200'
                return p1000_multi

            else:
                if current_tiprack != 'tips_1000':
                    protocol.move_labware(tips_1000, new_location='B3', use_gripper=True)
                    p1000_multi.configure_nozzle_layout(style=SINGLE, start="A1", tip_racks=[tips_1000])
                    current_tiprack = 'tips_1000'
                return p1000_multi

        # Perform sample and diluent transfers
        for i, row in unknown_samples.iterrows():
            source_well = sample_locations[i]
            destination_well = destination_wells[i]

            sample_vol = row['Sample Volume (µL)']
            diluent_vol = row['Diluent Volume (µL)']

            # Transfer sample
            if sample_vol > 0:
                pipette = use_pipette(sample_vol)
                pipette.pick_up_tip()
                pipette.transfer(
                    sample_vol,
                    temp_adapter[source_well],
                    plate3[destination_well],
                    rate=0.5,
                    new_tip='never'
                )
                pipette.drop_tip()

            # Transfer diluent
            if diluent_vol > 0:
                pipette = use_pipette(diluent_vol)
                pipette.pick_up_tip()
                pipette.transfer(
                    diluent_vol,
                    reservoir['A7'],
                    plate3[destination_well],
                    rate=0.5,
                    new_tip='never'
                )
                pipette.drop_tip()


    # ---------------- Click Reaction ----------------
    protocol.comment("Running click reaction")
    protocol.move_labware(labware=partial_50, new_location='B3', use_gripper=True)
    p50_multi.configure_nozzle_layout(style=SINGLE, start="A1", tip_racks=[partial_50]) #,
    
    #Pipette rhodamine azide (A3), tbta (A5), cuso4 (A2), and tcep (A4)
    p50_multi.transfer(1*(num_samples*2), 
                            temp_adapter['A3'], 
                            temp_adapter['A6'].bottom(z=0.1),
                            rate=speed,
                            mix_before=(1,10), 
                            delay=2,
                            disposal_vol=5,
                            #blow_out=True,
                            new_tip='always')

    p50_multi.transfer(3*(num_samples*2), 
                            temp_adapter['A5'], 
                            temp_adapter['A6'],
                            mix_before=(1,10),
                            rate=speed,
                            delay=3, 
                            new_tip='always')

    p50_multi.transfer(1*(num_samples*2), 
                            temp_adapter['A2'], 
                            temp_adapter['A6'], 
                            new_tip='always')

    p50_multi.transfer(1*(num_samples*2), 
                            temp_adapter['A4'], 
                            temp_adapter['A6'], 
                            mix_after=(3,30),
                            new_tip='always')

    # Make sure the click reagents are well mixed
    click_volume = 6*(final_volume/50)
    def mix_click_reagents():
        volume_click_reaction = final_volume + click_volume
        location = temp_adapter['A6']
        pipette = None
        moved_partial_50 = False
        moved_tips_1000 = False

        if volume_click_reaction < 100:
            positions_mixing = [1, 2, 3]
            pipette = p50_multi
        elif 100 < volume_click_reaction < 250:
            positions_mixing = [1, 4, 9]
            protocol.move_labware(labware=partial_50, new_location='B4', use_gripper=True)
            moved_partial_50 = True
            pipette = p1000_multi
        elif 250 < volume_click_reaction < 500:
            positions_mixing = [1, 6, 11]
            protocol.move_labware(labware=partial_50, new_location='B4', use_gripper=True)
            moved_partial_50 = True
            protocol.move_labware(labware=tips_1000, new_location='B3', use_gripper=True)
            moved_tips_1000 = True
            p1000_multi.configure_nozzle_layout(style=SINGLE, start="A1",tip_racks=[tips_1000])
            pipette = p1000_multi
        elif 500 < volume_click_reaction < 1000:
            positions_mixing = [1, 10, 16]
            protocol.move_labware(labware=partial_50, new_location='B4', use_gripper=True)
            moved_partial_50 = True
            protocol.move_labware(labware=tips_1000, new_location='B3', use_gripper=True)
            moved_tips_1000 = True
            p1000_multi.configure_nozzle_layout(style=SINGLE, start="A1",tip_racks=[tips_1000])
            pipette = p1000_multi
        else:
            positions_mixing = [1, 1, 1]  # Fallback
            pipette = p1000_multi

        # Perform mixing
        pipette.pick_up_tip()  # Optional if not done earlier
        pipette.aspirate(final_volume / 2, location.bottom(z=positions_mixing[0]))
        pipette.dispense(final_volume / 2, location.bottom(z=positions_mixing[1]))
        pipette.aspirate(final_volume / 3, location.bottom(z=positions_mixing[2]))
        pipette.dispense(final_volume / 3, location.bottom(z=positions_mixing[0]))
        pipette.mix(3, final_volume, location.bottom(z=positions_mixing[0]))
        pipette.drop_tip()

        # Move labware back if applicable
        if moved_tips_1000:
            protocol.move_labware(labware=tips_1000, new_location='C4', use_gripper=True)
        if moved_partial_50:
            protocol.move_labware(labware=partial_50, new_location='B3', use_gripper=True)

    # Call the function
    mix_click_reagents()


    # Pipette the click reaction premix
    protocol.move_labware(labware=partial_50, new_location='B3', use_gripper=True)
    click_volume = 6*(final_volume/50)
    p50_multi.transfer(click_volume, 
                            temp_adapter['A6'], 
                            [plate3[i] for i in destination_wells],
                            rate=speed-0.1,
                            delay=2,
                            mix_before=(3, 40),
                            mix_after=(3,30),
                            new_tip='always')

    # Step 11: shake the sample plate for click reaction
    protocol.move_labware(labware=plate3, new_location=heater_shaker, use_gripper=True)
    heater_shaker.close_labware_latch()
    heater_shaker.set_and_wait_for_shake_speed(1000)
    protocol.delay(minutes=90)
    heater_shaker.deactivate_shaker()
    heater_shaker.open_labware_latch()
    thermocycler.open_lid()
    protocol.move_labware(labware=plate3, new_location=thermocycler, use_gripper=True)

    # Add the loading buffer and move to the thermocylcer to seal and store.
    loading_buffer_volume = round((final_volume+click_volume) / 3, 1)
    columns = sorted(set(well[1:] for well in destination_wells), key=int)
    column_targets = [f'A{col}' for col in columns]
    p50_multi.configure_nozzle_layout(style=ALL, tip_racks=[partial_50])
    p50_multi.transfer(loading_buffer_volume, 
                            reservoir['A9'], 
                            [plate3[well] for well in column_targets], 
                            mix_after=(3, 40), 
                            new_tip='always')
    thermocycler.close_lid()
    thermocycler.set_block_temperature(95)  # Hold at 4°C
    protocol.delay(minutes=5)
    thermocycler.set_block_temperature(4)  # Hold at 4°C

