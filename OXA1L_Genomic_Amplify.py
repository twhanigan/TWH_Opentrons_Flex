from opentrons import protocol_api
from opentrons.protocol_api import SINGLE, ALL
import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
#import matplotlib.pyplot as plt
import subprocess
from pathlib import Path
import datetime

metadata = {
    'protocolName': 'BCA Assay with Normalization and Video Recording',
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
    
    num_samples = 7 #change this to the number of samples you need to run. The maximum is 18.
    # Change these if not using 96-well
    num_rows = 8  # A-H
    num_replicates = 3  # the number of replicates

    #Start recording the video
    video_process = subprocess.Popen(["python3", "/var/lib/jupyter/notebooks/record_video.py"])

    # Load modules
    heater_shaker = protocol.load_module('heaterShakerModuleV1', 'D1')
    thermocycler = protocol.load_module('thermocyclerModuleV2')
    temp_module = protocol.load_module('temperature module gen2', 'C1')
    mag_block = protocol.load_module('magneticBlockV1', 'D2')
    chute = protocol.load_waste_chute()

    # Load adapters
    temp_adapter = temp_module.load_labware('opentrons_24_aluminumblock_nest_1.5ml_screwcap')

    #set the temp module to 0c
    temp_module.set_temperature(celsius=10)
    
    # Load labware
    tips_50 = protocol.load_labware('opentrons_flex_96_filtertiprack_50ul', 'A4')
    partial_50 = protocol.load_labware(load_name="opentrons_flex_96_filtertiprack_50ul",location="A3")
    tips_200 = protocol.load_labware('opentrons_flex_96_filtertiprack_200ul', 'B4')
    partial_200 = protocol.load_labware(load_name="opentrons_flex_96_filtertiprack_200ul",location="B3")
    reservoir = protocol.load_labware('nest_12_reservoir_15ml', 'C2')

    # Load PCR plate on thermocycler
    pcr_plate = thermocycler.load_labware('opentrons_96_wellplate_200ul_pcr_full_skirt')
    
    # Load reagent tubes on heater shaker with adapter
    temp_module = heater_shaker.load_labware(
        'opentrons_24_aluminumblock_nest_2ml_screwcap')
    
    # Liquid definitionss
    OXA1L_F_Primer = protocol.define_liquid(name = 'OXA1L_F_Primer', display_color="#704848",)
    OXA1L_R_Primer = protocol.define_liquid(name = 'OXA1L_R_Primer', display_color="#704300",)
    HF_Buffer = protocol.define_liquid(name = 'HF_Buffer', display_color="#704900",)
    Phusion = protocol.define_liquid(name = 'Phusion', display_color="#701100",)
    ddH2O = protocol.define_liquid(name = 'ddH2O', display_color="#704900",)
    DMSO = protocol.define_liquid(name = 'DMSO', display_color="#704300",)
    dNTPs = protocol.define_liquid(name='dNTPs', display_color="#FFC0CB")  # Pink

    sample_liquids = [protocol.define_liquid(name = f'Sample {i + 1}', display_color="#FFA000",) for i in range(num_samples)]

    temp_adapter['A1'].load_liquid(liquid=OXA1L_F_Primer, volume=1500) #click
    temp_adapter['A2'].load_liquid(liquid=OXA1L_R_Primer, volume=1500) #click
    temp_adapter['A3'].load_liquid(liquid=HF_Buffer, volume=1500) #click
    temp_adapter['A4'].load_liquid(liquid=Phusion, volume=1500) #click
    temp_adapter['A5'].load_liquid(liquid=ddH2O, volume=1500) #click
    temp_adapter['A6'].load_liquid(liquid=DMSO, volume=1500) #click
    temp_adapter['A7'].load_liquid(liquid=dNTPs, volume=1500) #click

    # Load pipettes
    p50_multi = protocol.load_instrument('flex_8channel_50', 'left', tip_racks=[tiprack_50])
    p1000_multi = protocol.load_instrument('flex_8channel_1000', 'right', tip_racks=[tiprack_200])
    
    # Open thermocycler lid
    thermocycler.open_lid()
    
    # Transfer reagents
    # 1. Transfer OXA1L_F (2.5 μL)
    p50_multi.transfer(2.5, temp_module['A1'], pcr_plate.columns()[0], new_tip='once')
    
    # 2. Transfer OXA1L_R (2.5 μL)
    p50_multi.transfer(2.5, temp_module['A2'], pcr_plate.columns()[0], new_tip='once')
    
    # 3. Transfer Phusion buffer (10 μL)
    p50_multi.transfer(10, temp_module['A3'], pcr_plate.columns()[0], new_tip='once')
    
    # 4. Transfer dNTP mix (1 μL)
    p50_multi.transfer(1, temp_module['A4'], pcr_plate.columns()[0], new_tip='once')
    
    # 5. Transfer Phusion polymerase (0.5 μL)
    p50_multi.transfer(0.5, temp_module['A5'], pcr_plate.columns()[0], new_tip='once')
    
    # 6. Transfer H2O (31.5 μL)
    p50_multi.transfer(31.5, temp_module['A6'], pcr_plate.columns()[0], new_tip='once')

    # Close thermocycler lid and set temperature
    thermocycler.close_lid()
    
    # 7. Run PCR cycles
    # Initial denaturation
    thermocycler.execute_profile(
        steps=[{'temperature': 98, 'hold_time_seconds': 15}],
        repetitions=1,
        block_max_volume=50
    )
    
    # 34 cycles of PCR
    thermocycler.execute_profile(
        steps=[
            {'temperature': 98, 'hold_time_seconds': 10},
            {'temperature': 59, 'hold_time_seconds': 15},
            {'temperature': 72, 'hold_time_seconds': 150}
        ],
        repetitions=34,
        block_max_volume=50
    )
    
    # Final extension
    thermocycler.execute_profile(
        steps=[{'temperature': 72, 'hold_time_seconds': 600}],
        repetitions=1,
        block_max_volume=50
    )
    
    # Hold at 4°C
    thermocycler.set_block_temperature(4)
    
    # Open lid
    thermocycler.open_lid()