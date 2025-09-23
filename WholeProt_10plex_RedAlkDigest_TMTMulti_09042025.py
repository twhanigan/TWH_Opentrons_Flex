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
    'protocolName': 'TMT-labeled whole proteome sample preparation (No BCA) 09092025',
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
    video_process = subprocess.Popen(["python3", "/var/lib/jupyter/notebooks/record_video_wholeprot.py"])

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
    tips_200 = protocol.load_labware(load_name='opentrons_flex_96_filtertiprack_200ul', location='B4')
    tips_1000 = protocol.load_labware(load_name="opentrons_flex_96_filtertiprack_1000ul",location="B3")
    tips_1000_2 = protocol.load_labware(load_name="opentrons_flex_96_filtertiprack_1000ul",location="D4")
    plate3 = protocol.load_labware('thermoscientificnunc_96_wellplate_2000ul', location="B2") 
    plate4 = thermocycler.load_labware('opentrons_96_wellplate_200ul_pcr_full_skirt')
    tmt_plate = protocol.load_labware('opentrons_96_wellplate_200ul_pcr_full_skirt', location='A2')
    reservoir = protocol.load_labware('nest_12_reservoir_15ml', 'C2')
    tiprack_locations = {"tips_50": "A3", "tips_200": "B4", "tips_1000": "B3", "tips_1000_2": "D4"}

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
    ethanol_80 = protocol.define_liquid(name='80% Ethanol', display_color="#4682B4")  # SteelBlue
    trypsin = protocol.define_liquid(name='Trypsin in EPPS', display_color="#9900FF")  # Gold
    cacl2 = protocol.define_liquid(name='CaCl2', display_color="#FF3300")  # LimeGreen
    hydroxylamine = protocol.define_liquid(name='hydroxylamine', display_color="#8A2BE2")   # Blue Violet
    sample_liquids = [protocol.define_liquid(name = f'Sample {i + 1}', display_color="#FFA000",) for i in range(protocol.params.num_samples)]

    # Reservoir assignments for washes and digestion
    reservoir['A1'].load_liquid(liquid=bsa_reag_a, volume=15000)  
    reservoir['A2'].load_liquid(liquid=bsa_reag_b, volume=15000)  
    reservoir['A3'].load_liquid(liquid=bsa_reag_c, volume=15000)  
    reservoir['A4'].load_liquid(liquid=excess_lysis, volume=15000) 
    reservoir['A7'].load_liquid(liquid=abs_ethanol, volume=15000)
    reservoir['A8'].load_liquid(liquid=ethanol_80,volume=15000)
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

    # ---------------- Normalizing BCA Assay ----------------
    protocol.comment("Place BCA assay absorbance data in /var/lib/jupyter/notebooks/Data, load new deep well plate into hs_adapter, and new tube rack into A2 (with click mix, beads and TCEP/IAA)")
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

    #Define the directory path
    directory = Path("/var/lib/jupyter/notebooks/Data/")

    # Get today's date in YYMMDD format and find file with date
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

    # Create a list of well names (A1 to H12) and flatten to a dataframe
    well_names = [f"{row}{col}" for col in range(1, 13) for row in "ABCDEFGH"]
    absorbance_values = df.values.flatten()
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
    
    # Write the output and image of data plot to the instrument jupyter notebook directory
    normalized_samples = unknown_samples[['Sample', 'Protein Concentration (mg/mL)', 'Sample Volume (µL)', 'Diluent Volume (µL)']].reset_index().drop(columns='index')
    filename = f"Protocol_output_{today_date}.csv"
    output_file_destination_path = directory.joinpath(filename)
    normalized_samples.to_csv(output_file_destination_path)
    rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    destination_wells  = [f'{rows[i % 8]}{(i // 8)+ 1}' for i in range(len(normalized_samples))]

    # Dilute sample in lysis buffer to 1 mg/ml on deep well plate
    destination_wells  = [f'{rows[i % 8]}{(i // 8)+ 1}' for i in range(len(normalized_samples))]
    tiprack_locations = {"tips_50": "A3", "tips_200": "B4", "tips_1000": "B3", "tips_1000_2": "D4"}

    def normalize_samples():
        fv = protocol.params.final_volume
        pipette = None
        # ≤50 µL  → p50 + 50 µL tips
        if fv <= 50:
            if tiprack_locations["tips_1000"] != "C4":
                protocol.move_labware(tips_1000, new_location="C4", use_gripper=True)
                tiprack_locations["tips_1000"] = "C4"
            p50_multi.configure_nozzle_layout(style=SINGLE, start="A1", tip_racks=[tips_50])
            pipette = p50_multi
        # 50–200 µL → p1000 + 200 µL tips  (per your request)
        elif 50 < fv <= 200:
            if tiprack_locations["tips_1000"] != "C4":
                protocol.move_labware(tips_1000, new_location="C4", use_gripper=True)
                tiprack_locations["tips_1000"] = "C4"
            if tiprack_locations["tips_200"] != "B3":
                protocol.move_labware(tips_200, new_location="B3", use_gripper=True)
                tiprack_locations["tips_200"] = "B3"
            p1000_multi.configure_nozzle_layout(style=SINGLE, start="A1", tip_racks=[tips_200])
            pipette = p1000_multi
        # >200 µL → p1000 + 1000 µL tips
        else:
            if tiprack_locations["tips_1000"] != "B3":
                protocol.move_labware(tips_1000, new_location="B3", use_gripper=True)
                tiprack_locations["tips_1000"] = "B3"
            p1000_multi.configure_nozzle_layout(style=SINGLE, start="A1", tip_racks=[tips_1000])
            pipette = p1000_multi

        # one tip per sample: sample + diluent to same destination
        for i, row in normalized_samples.iterrows():
            source_well = sample_locations[i]
            destination_well = destination_wells[i]
            normalized_volume = float(row['Sample Volume (µL)'])
            diluent_volume = float(fv - normalized_volume)

            pipette.transfer(
                normalized_volume,
                temp_adapter[source_well],
                plate3[destination_well].bottom(z=0.1),
                rate=speed,
                new_tip='always'
            )
            pipette.transfer(
                diluent_volume,
                reservoir['A4'],   # diluent
                plate3[destination_well].bottom(z=1),
                rate=speed,
                new_tip='always',
                mix_after=(3, 10)
            )

    normalize_samples()

    # ---------------- Reduction, Alkylation and Digestion ----------------
    #set the heater_shaker temp to 37 c for reduce
    heater_shaker.set_and_wait_for_temperature(37)

    #Configure the p1000 pipette to use single tip NOTE: this resets the pipettes tip racks!
    if tiprack_locations["tips_200"] != "B4":
        protocol.move_labware(tips_200, new_location="B4", use_gripper=True)
        tiprack_locations["tips_200"] = "B4"
    if tiprack_locations["tips_1000"] != "C4":
        protocol.move_labware(tips_1000, new_location="C4", use_gripper=True)
        tiprack_locations["tips_1000"] = "C4"

    #Pipette mix K2CO3, and TCEP and transfer to samples above the 200 uL liquid miniscus.
    Redu_mix = protocol.params.num_samples*5*2/2 # 10 mM TCEP, 30 mM K2CO2 final conc from a 200 mM/600 mM stock
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
    heater_shaker.deactivate_heater()
    protocol.comment('Waiting for heater_shaker to cool')
    protocol.delay(minutes=3)
    heater_shaker.open_labware_latch()

    # add IAA and alkylate for 30 minutes
    protocol.move_labware(labware=plate3, new_location='B2', use_gripper=True)
    p50_multi.distribute(7, 
                            temp_adapter['A5'], 
                            [plate3[i].bottom(z=1) for i in destination_wells],
                            rate=speed,
                            mix_before=(1,50),
                            new_tip='always')

    protocol.move_labware(labware=plate3, new_location=heater_shaker, use_gripper=True)
    heater_shaker.close_labware_latch()
    heater_shaker.set_and_wait_for_shake_speed(1000)
    protocol.delay(minutes=30)
    heater_shaker.deactivate_shaker()
    heater_shaker.open_labware_latch()
    protocol.move_labware(labware=plate3, new_location='B2', use_gripper=True)

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
                            rate=speed,
                            delay=2)
    
    # incubate beads with shaking for 5 minutes
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
                            reservoir['A7'], 
                            [well.top(z=-5) for well in destination_wells_col], 
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
        p1000_multi.aspirate(200, well.bottom(z=2), rate=0.1)
        protocol.delay(seconds=2)
        p1000_multi.blow_out(chute)
    p1000_multi.drop_tip()  

    # Remove remainder of EtOH and Add 80% EtOH and wash the beads three times with 200 uL 80% EtOH
    p50_multi.configure_nozzle_layout(style=ALL, tip_racks=[tips_50])
    for i in range(3):
        # Ensure p1000 has a tip before any aspirate/distribute with new_tip='never'
        if not p1000_multi.has_tip:
            p1000_multi.pick_up_tip()

        # --- REMOVE SUPERNATANT (slow, above pellet) ---
        for well in destination_wells_col:
            p1000_multi.aspirate(180, well.bottom(z=1), rate=0.1)
            protocol.delay(seconds=2)
            p1000_multi.blow_out(chute)

        # --- ADD 80% EtOH (uses same tip; new_tip='never' assumes tip is on) ---
        p1000_multi.distribute(
            180,
            reservoir['A8'],
            [w.top(z=-10) for w in destination_wells_col],
            rate=speed,
            new_tip='never'
        )

        # Shake, then back to magnet
        protocol.move_labware(labware=plate3, new_location=heater_shaker, use_gripper=True)
        heater_shaker.close_labware_latch()
        heater_shaker.set_and_wait_for_shake_speed(1000)
        protocol.delay(seconds=60)
        heater_shaker.deactivate_shaker()
        heater_shaker.open_labware_latch()
        protocol.move_labware(labware=plate3, new_location=mag_block, use_gripper=True)
        protocol.delay(minutes=3)

        # After the final (3rd) wash, pull the last ~180 and residual with p50/account for residual
        if i == 2:
            if not p1000_multi.has_tip:
                p1000_multi.pick_up_tip()
            for well in destination_wells_col:
                p1000_multi.aspirate(180, well.bottom(z=1), rate=0.1)
                protocol.delay(seconds=2)
                p1000_multi.blow_out(chute)
            if p1000_multi.has_tip:
                p1000_multi.drop_tip()
            if not p50_multi.has_tip:
                p50_multi.pick_up_tip()
            for well in destination_wells_col:
                p50_multi.aspirate(50, well.bottom(z=0.1), rate=0.1)
                protocol.delay(seconds=2)
                p50_multi.blow_out(chute)
            p50_multi.drop_tip()

    # If you kept the same p1000 tip through all washes, drop it at the end:
    if p1000_multi.has_tip:
        p1000_multi.drop_tip()

    # Air dry
    protocol.delay(minutes=5)

    # Resuspend in 2 M urea in EPPS for digestion setup
    protocol.move_labware(labware=plate3, new_location='B2', use_gripper=True)
    p1000_multi.transfer(
        47.5,
        reservoir['A12'],
        [well for well in destination_wells_col],
        mix_before=(1, 50),
        mix_after=(3, 30),
        new_tip='always'
    )

    # ---------------- Digestion ----------------
    #Configure the p1000 pipette to use All tips 
    protocol.move_labware(labware=tips_1000_2, new_location="D4", use_gripper=True)
    p1000_multi.configure_nozzle_layout(style=SINGLE, start='A1' ,tip_racks=[tips_1000])

    # Add CaCl2, trypsin in epps, and move to shaker
    p1000_multi.distribute(2.5, 
                            temp_adapter['D5'], 
                            [plate3[i].bottom(z=0.5) for i in destination_wells],
                            disposal_vol=10, 
                            new_tip='once')

    # Add Trypsin. For 50ug protein, resuspend trypsin in 2000 uL 2 M urea EPPS and place 500 uL on the robot. 
    p1000_multi.distribute(50, 
                            temp_adapter['D6'], 
                            [plate3[i].top(z=-10) for i in destination_wells],
                            disposal_vol=0,
                            rate=speed,
                            mix_before=(1,150),
                            #mix_after=(3,500),
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

    # Remove samples from beads on plate3
    protocol.move_labware(labware=tips_200, new_location="C3", use_gripper=True)
    p1000_multi.configure_nozzle_layout(style=ALL, tip_racks=[tips_200])
    rows = "ABCDEFGH"  # ensure defined
    src_cols_count = (protocol.params.num_samples + 7) // 8  
    dest_start_col_idx = src_cols_count  # next free column to the right (0-based)
    max_cols = 12
    if dest_start_col_idx + src_cols_count > max_cols:
        raise ValueError("Not enough columns on plate3 to copy samples to the right.")

    for c in range(src_cols_count):
        src_top = plate3.columns()[c][0]                               
        dst_top = plate3.columns()[dest_start_col_idx + c][0]

        p1000_multi.transfer(
            100,                               
            src_top.bottom(z=0.1),
            dst_top.bottom(z=0.1),
            new_tip='always',
            rate=0.2,
        )

    # ---------------- Post-Digestion TMT Reaction ----------------
    protocol.comment("Running TMT reaction")
    total_protein = protocol.params.final_volume * protocol.params.target_concentration  # µg
    tmt_vol = 10 * (47.5 + 2.5 + 50) / total_protein  # µL for 10 µg

    # Move 10 ug of sample to plate4 from the plate3 wells containing sample without beads
    protocol.move_labware(labware=tips_200, new_location="B4", use_gripper=True)
    protocol.move_labware(labware=plate3, new_location='B2', use_gripper=True)

    p50_multi.configure_nozzle_layout(style=ALL, tip_racks=[tips_50])

    for c in range(src_cols_count):  # includes the final partial column
        src_anchor = plate3.columns()[dest_start_col_idx + c][0]  # A{src_col}
        dst_anchor = plate4.columns()[c][0]                       # A{dst_col}
        p50_multi.transfer(
            tmt_vol,
            src_anchor.bottom(z=0.5),
            dst_anchor.bottom(z=0.1),
            new_tip='always',
            rate=speed,
        )
   

    # Add TMT reagent from the row-wise tmt_plate sources to matching plate4 wells
    p50_multi.configure_nozzle_layout(style=SINGLE, start="A1", tip_racks=[tips_50])
    tmt_dilution_wells = [plate4.wells_by_name()[f"{rows[i % 8]}{(i // 8) + 1}"] for i in range(protocol.params.num_samples)]
    tmt_sources = [tmt_plate.wells_by_name()[f"{protocol.params.tmt_row}{c}"] for c in range(1, protocol.params.num_samples + 1)]

    for source, dest in zip(tmt_sources, tmt_dilution_wells):
        p50_multi.transfer(
            5,
            source,
            dest.bottom(z=0.1),
            disposal_volume=0,
            rate=0.2,
            new_tip='always',
            mix_after=(3, 10)
        )

    # React
    thermocycler.close_lid()
    protocol.delay(minutes=120)
    thermocycler.open_lid()

    # Quench with hydroxylamine
    for dest in tmt_dilution_wells:
        p50_multi.transfer(
            1,
            temp_adapter['D4'],
            dest.bottom(z=0.1),
            rate=speed,
            mix_before=(1,5),
            mix_after=(2, 5),
            new_tip='always'   # avoid cross-contamination
        )

    thermocycler.close_lid()
    protocol.delay(minutes=5)
    video_process.terminate()
    thermocycler.set_block_temperature(4)  # hold at 4°C
