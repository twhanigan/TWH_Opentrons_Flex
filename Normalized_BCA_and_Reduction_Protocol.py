from opentrons import protocol_api
from opentrons.protocol_api import SINGLE, ALL
import pandas as pd
import numpy as np
import subprocess
from pathlib import Path
import datetime

metadata = {
    'protocolName': 'Normalized BCA and Reduction Workflow',
    'author': 'Assistant Rewrite',
    'description': 'Starts at BCA normalization using relocated labware.'
}

requirements = {
    "robotType": "Flex",
    "apiLevel": "2.21"
}

def run(protocol: protocol_api.ProtocolContext):
    protocol.comment("Start: Normalizing BCA Assay")

    # Load modules
    heater_shaker = protocol.load_module('heaterShakerModuleV1', 'D1')
    mag_block = protocol.load_module('magneticBlockV1', 'D2')
    temp_module = protocol.load_module('temperature module gen2', 'C1')
    chute = protocol.load_waste_chute()

    # Load adapters
    plate_adapter = heater_shaker.load_adapter('opentrons_96_deep_well_adapter')
    temp_adapter = temp_module.load_labware('opentrons_24_aluminumblock_nest_1.5ml_screwcap')

    # Load labware (with updated locations from previous moves)
    plate3 = protocol.load_labware('nest_96_wellplate_2ml_deep', 'C3')
    epp_rack = protocol.load_labware('opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap', 'B2')
    reservoir = protocol.load_labware('nest_12_reservoir_15ml', 'C2')

    tips_50_partial = protocol.load_labware('opentrons_flex_96_filtertiprack_50ul', 'A2')
    tips_1000_partial = protocol.load_labware('opentrons_flex_96_filtertiprack_1000ul', 'B3')
    tips_1000_complete = protocol.load_labware('opentrons_flex_96_filtertiprack_1000ul', 'C4')
    tips_1000_2 = protocol.load_labware('opentrons_flex_96_filtertiprack_1000ul', 'A4')
    tips_200 = protocol.load_labware('opentrons_flex_96_filtertiprack_200ul', 'B3')

    # Load instruments with updated tip racks
    p50_multi = protocol.load_instrument('flex_8channel_50', 'left')
    p1000_multi = protocol.load_instrument('flex_8channel_1000', 'right')

    # Configure single channel usage
    p50_multi.configure_nozzle_layout(style=SINGLE, start="A1", tip_racks=[tips_50_partial])
    p1000_multi.configure_nozzle_layout(style=SINGLE, start="A1", tip_racks=[tips_1000_partial])

    # Pause for user to upload file
    protocol.pause("Place BCA assay absorbance data in /var/lib/jupyter/notebooks/TWH and press resume")

    # Locate file using external script
    find_file = subprocess.Popen(
        ['python3', "/var/lib/jupyter/notebooks/wait_for_file.py"],
        stdout=subprocess.PIPE, text=True
    )
    stdout, stderr = find_file.communicate()

    if stderr:
        raise ValueError(f"Error while waiting for file: {stderr}")
    file_path = stdout.splitlines()[1]
    if not file_path:
        raise ValueError("No file path returned by wait_for_file.py")

    protocol.comment(f"Successfully loaded: {file_path}")

    # Read and process absorbance data
    df = pd.read_excel(file_path, header=5, nrows=8, usecols="C:N")
    well_names = [f"{row}{col}" for col in range(1, 13) for row in "ABCDEFGH"]
    absorbance_values = df.values.flatten()
    initial_df = pd.DataFrame({'Well': well_names, 'Absorbance': absorbance_values})

    # Normalize data
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

    samples_1_to_8 = final_df.iloc[:8]
    samples_1_to_8['Mean Absorbance'] = samples_1_to_8[['Replicate 1', 'Replicate 2', 'Replicate 3']].mean(axis=1)
    protein_concentrations = [10, 5, 2.5, 1.25, 0.625, 0.3125, 0.15625, 0]
    samples_1_to_8['Protein Concentration (mg/mL)'] = protein_concentrations

    slope, intercept = np.polyfit(samples_1_to_8['Protein Concentration (mg/mL)'], samples_1_to_8['Mean Absorbance'], 1)
    unknown_samples = final_df.iloc[8:]
    unknown_samples['Mean Absorbance'] = unknown_samples[['Replicate 1', 'Replicate 2', 'Replicate 3']].mean(axis=1)
    unknown_samples['Protein Concentration (mg/mL)'] = (unknown_samples['Mean Absorbance'] - intercept) / slope

    final_volume = 0.440  # mL
    target_concentration = 1  # mg/mL
    unknown_samples['Sample Volume (mL)'] = (target_concentration * final_volume) / unknown_samples['Protein Concentration (mg/mL)']
    unknown_samples['Diluent Volume (mL)'] = final_volume - unknown_samples['Sample Volume (mL)']
    unknown_samples.loc[unknown_samples['Sample Volume (mL)'] > final_volume, ['Sample Volume (mL)', 'Diluent Volume (mL)']] = [final_volume, 0]

    today_date = datetime.date.today().strftime("%y%m%d")
    output_path = Path("/var/lib/jupyter/notebooks/TWH/") / f"Protocol_output_{today_date}.csv"
    unknown_samples.to_csv(output_path)

    protocol.comment("Normalized unknowns to 1 mg/mL.")

    # Dilute samples on deepwell plate
    num_samples = len(unknown_samples)
    rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    destination_wells = [f'{rows[i % 8]}{(i // 8) + 1}' for i in range(num_samples)]
    sample_locations = [f'B{i+1}' if i < 6 else f'C{i-5}' if i < 12 else f'D{i-11}' for i in range(num_samples)]

    for i, row in unknown_samples.iterrows():
        source_well = sample_locations[i]
        normalized_volume = row['Sample Volume (mL)']
        diluent_volume = 500 - normalized_volume
        dest = destination_wells[i]
        p1000_multi.transfer(diluent_volume, reservoir['A7'], plate3[dest], rate=0.5, new_tip='once')
        p50_multi.transfer(normalized_volume, temp_adapter[source_well], plate3[dest], rate=0.5, new_tip='once')

    # Desalting, reduction, and digestion
    p1000_multi.pick_up_tip()
    p1000_multi.mix(4, 200, epp_rack['B1'])
    p1000_multi.drop_tip()
    p1000_multi.distribute(60, epp_rack['B1'], [plate3[i] for i in destination_wells])
    added_vol = 60
    heater_shaker.close_labware_latch()
    heater_shaker.set_and_wait_for_shake_speed(1000)
    protocol.delay(minutes=5)
    protocol.move_labware(labware=tips_1000_partial, new_location="B4", use_gripper=True)
    protocol.move_labware(labware=tips_1000_complete, new_location="B3", use_gripper=True)
    p1000_multi.configure_nozzle_layout(style=ALL, tip_racks=[tips_1000_complete])
    num_full_columns = (num_samples + 7) // 8
    destination_columns = plate3.columns()[:num_full_columns]
    destination_wells_col = [col[0] for col in destination_columns]
    p1000_multi.distribute(600, reservoir['A9'], destination_wells_col, new_tip='once', air_gap=10)
    added_vol += 600
    heater_shaker.set_and_wait_for_shake_speed(1000)
    protocol.delay(minutes=5)
    heater_shaker.deactivate_shaker()
    heater_shaker.open_labware_latch()
    protocol.move_labware(labware=plate3, new_location=mag_block, use_gripper=True)
    protocol.delay(minutes=3)
    p1000_multi.consolidate(900, [well.bottom(z=0.2) for well in destination_wells_col], reservoir['A11'], rate=0.5, new_tip='once')
    added_vol -= 900

    # Reduction and alkylation
    protocol.move_labware(labware=plate3, new_location=plate_adapter, use_gripper=True)
    heater_shaker.close_labware_latch()
    p1000_multi.distribute(200, reservoir['A12'], destination_wells_col, mix_after=(3, 150), new_tip='once')
    added_vol += 200
    heater_shaker.set_and_wait_for_temperature(37)
    protocol.move_labware(labware=tips_1000_complete, new_location="C4", use_gripper=True)
    protocol.move_labware(labware=tips_1000_partial, new_location="B3", use_gripper=True)
    p1000_multi.configure_nozzle_layout(style=SINGLE, start="A1", tip_racks=[tips_1000_partial])
    redu_mix = num_samples * 50 * 1.2 / 2
    p1000_multi.transfer(redu_mix, epp_rack['C1'], epp_rack['C4'], new_tip='always')
    p1000_multi.transfer(redu_mix, epp_rack['C2'], epp_rack['C4'], new_tip='always')
    p1000_multi.distribute(50, epp_rack['C4'], [plate3[i] for i in destination_wells], new_tip='always')
    added_vol += 50
    heater_shaker.set_and_wait_for_shake_speed(1000)
    protocol.delay(minutes=30)
    heater_shaker.deactivate_shaker()
    p1000_multi.distribute(70, epp_rack['C3'], [plate3[i] for i in destination_wells], new_tip='always')
    heater_shaker.set_and_wait_for_shake_speed(1000)
    protocol.delay(minutes=30)
    heater_shaker.deactivate_shaker()
    added_vol += 70

    # Trypsin digestion
    protocol.move_labware(labware=tips_1000_2, new_location="C3", use_gripper=True)
    p1000_multi.configure_nozzle_layout(style=ALL, tip_racks=[tips_1000_2])
    p1000_multi.distribute(2.5, temp_adapter['D5'], [plate3[i] for i in destination_wells], new_tip='once')
    p1000_multi.distribute(150, temp_adapter['D6'], [plate3[i] for i in destination_wells], new_tip='once')
    added_vol += 2.5 + 150
    heater_shaker.close_labware_latch()
    heater_shaker.set_and_wait_for_shake_speed(1000)
    protocol.delay(minutes=960)
    heater_shaker.deactivate_shaker()
    heater_shaker.open_labware_latch()
    protocol.comment("Trypsin digestion complete. Ready for downstream processing.")
