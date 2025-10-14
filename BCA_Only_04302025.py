from opentrons import protocol_api
from opentrons.protocol_api import SINGLE, ALL
import pandas as pd
import numpy as np
import subprocess
from pathlib import Path
#import matplotlib.pyplot as plt
import datetime
import time

metadata = {
    'protocolName': 'BCA Assay w/out Normalization',
    'author': 'Assistant',
    'description': 'Serial dilution of BSA standard and sample processing. This includes cooling samples to 4c, heating plate to 37c with shaking and recording a video of the whole process. '
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


def run(protocol: protocol_api.ProtocolContext):
    protocol.comment(
        "Place BSA Standard in A1, reagents a, b, c and Lysis buffer in reservoir, and samples starting in row B")
    protocol.comment("Running the BCA assay")

    #Speed of pipetting NP40 lysis buffer=0.35, 2M Urea in EPPS=0.75
    speed= 0.3
    num_rows = 8  # A-H

    #Start recording the video
    video_process = subprocess.Popen(["python3", "/var/lib/jupyter/notebooks/record_video_bca.py"])

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
    heater_shaker.set_and_wait_for_temperature(50)

    #set the temp module to 0c
    temp_module.set_temperature(celsius=10)
    
    # Load labware
    partial_50 = protocol.load_labware(load_name="opentrons_flex_96_filtertiprack_50ul",location="B3")
    tips_200 = protocol.load_labware(load_name="opentrons_flex_96_filtertiprack_200ul",location="A3")
    tips_1000 = protocol.load_labware('opentrons_flex_96_filtertiprack_1000ul', 'C4')
    plate1 = protocol.load_labware('opentrons_96_wellplate_200ul_pcr_full_skirt', 'A2')
    plate2 = protocol.load_labware('corning_96_wellplate_360ul_flat', 'B2')
    reservoir = protocol.load_labware('nest_12_reservoir_15ml', 'C2')
    
    # Liquid definitions
    bsa_standard = protocol.define_liquid(name = 'BSA Standard', display_color="#704848",)
    bsa_reag_a = protocol.define_liquid(name = 'Reagent A', display_color="#704300",)
    bsa_reag_b = protocol.define_liquid(name = 'Reagent B', display_color="#704900",)
    bsa_reag_c = protocol.define_liquid(name = 'Reagent C', display_color="#701100",)
    excess_lysis = protocol.define_liquid(name='Excess Lysis Buffer', display_color="#FFC0CB")  # Pink
    sample_liquids = [protocol.define_liquid(name = f'Sample {i + 1}', display_color="#FFA000",) for i in range(protocol.params.num_samples)]

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

    # Steps 1: Add lysis buffer to column 1 of plate1. 
    p1000_multi.distribute(50, 
         reservoir['A7'],
         plate1[f'A{protocol.params.standards_col}'],
         rate = speed,
         delay = 2,
         new_tip='once')

    #Step 3: Configure the p50 pipette to use single tip NOTE: this resets the pipettes tip racks!
    p50_multi.configure_nozzle_layout(style=SINGLE, start="A1",tip_racks=[partial_50])

    # Step 4: Transfer BSA standard (20 mg/ml) to first well of column 1
    p50_multi.transfer(50,
        temp_adapter['A1'],
        plate1[f'A{protocol.params.standards_col}'],
        rate = speed,
        delay = 2,
        mix_after=(3, 40),
        new_tip='once',
        disposal_vol=0)

    # Step 5: Perform serial dilution down column 1
    rows = ['A','B', 'C', 'D', 'E', 'F', 'G']
    p50_multi.pick_up_tip()
    for source, dest in zip(rows[:-1], rows[1:]):
        p50_multi.transfer(50,
                         plate1[f'{source}{protocol.params.standards_col}'],
                         plate1[f'{dest}{protocol.params.standards_col}'],
                         rate = speed,
                         mix_after=(3, 40),
                         new_tip='never', 
                         disposal_vol=0)
    p50_multi.drop_tip()

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
    p50_multi.configure_nozzle_layout(style=ALL, tip_racks=[partial_50]) #, 

    #Step 10: Pipette triplicate of controls from plate1 column 1 to plate2 columns 1,2,3 
    p50_multi.distribute(5, 
                        plate1[f'A{protocol.params.standards_col}'], 
                        [plate2[f'A{i}'].bottom(z=0.1) for i in range(1, 4)],
                        rate= speed,
                        disposal_vol=5)

    # Step 11: move the 50 uL partial tips to C3 and the 200uL complete tips to B3
    protocol.move_labware(labware=tips_1000, new_location="C3", use_gripper=True)

    #Step 12: Load the p1000 with full tip rack
    p1000_multi.configure_nozzle_layout(style=ALL, tip_racks=[tips_1000]) #,

    # Step 13: Add reagent A
    p1000_multi.distribute(50,
                        reservoir['A1'],
                        plate2.wells(),
                        new_tip='once',
                        disposal_vol=50,)

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
                        disposal_vol=5)

    #Step 16: move plate 2 to the heater shaker and incubate at 37c
    heater_shaker.open_labware_latch()
    protocol.move_labware(labware=plate2, new_location=hs_adapter,use_gripper=True)
    heater_shaker.close_labware_latch()
    heater_shaker.set_and_wait_for_shake_speed(500)
    protocol.delay(minutes=20)

    #Step 17 deactivate heater shaker and temp modules
    heater_shaker.deactivate_shaker()
    heater_shaker.deactivate_heater()
    heater_shaker.open_labware_latch()
    # Stop video recording after the main task is completed
    video_process.terminate()
