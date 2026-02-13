from opentrons import protocol_api
from opentrons.protocol_api.labware import Labware
from opentrons.protocol_api import SINGLE, ALL
import pandas as pd
import numpy as np
import subprocess
from pathlib import Path
import math
import datetime
import time

metadata = {
    'protocolName': 'TMT-Based Photolabeling and Chemical Proteomics Dummy Data',
    'author': 'Om Patel and Thomas Hanigan',
    'description': 'Sample normalization, click reaction, and denaturation for Gel-based photolabeling.'
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
        default=440,
        minimum=100,
        maximum=1000,
        unit="µL"
    )

def run(protocol: protocol_api.ProtocolContext):
    # ---------------- BCA Assay for Preprepared Lysates ----------------
    protocol.comment("Running the BCA assay")
    
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
    plate3 = protocol.load_labware('thermoscientificnunc_96_wellplate_2000ul', location='A4')  # original location (D3 reserved for waste chute) 
    plate4 = protocol.load_labware('opentrons_96_wellplate_200ul_pcr_full_skirt', location=protocol_api.OFF_DECK)
    tmt_plate = protocol.load_labware('opentrons_96_wellplate_200ul_pcr_full_skirt', location=protocol_api.OFF_DECK)
    reservoir = protocol.load_labware('nest_12_reservoir_15ml', 'C2')
    tiprack_locations = {
        "tips_50": "A3",
        "tips_200": "B3",
        "tips_1000": "C4",
        "tips_1000_2": "D4"
    }

    # --- Helper: standardize column-block addressing on plate3 ---
    # Samples are assumed to be packed down columns starting at column 1 (A1–H1, then A2–H2, ...).
    # block_index=0 -> columns 1..cols_used; block_index=1 -> next cols_used columns; etc.
    def get_plate3_block(block_index: int):
        cols_used = math.ceil(protocol.params.num_samples / num_rows)
        start = block_index * cols_used
        end = start + cols_used

        if end > len(plate3.columns()):
            raise ValueError(
                f"Not enough columns on plate3 for block_index={block_index}. "
                f"cols_used={cols_used}, need up to column {end}, plate has {len(plate3.columns())} columns."
            )

        cols = plate3.columns()[start:end]
        wells_col = [c[0] for c in cols]  # A-row anchors for column-based transfers
        wells_single = [w for c in cols for w in c][:protocol.params.num_samples]  # A->H per col, truncated
        return cols_used, cols, wells_col, wells_single
    
    # load the liquids to the tube racks and reservoirs (all unique + vivid)
    bsa_standard=protocol.define_liquid(name='BSA Standard',display_color="#8B4513")  # brown (saddlebrown)
    bsa_reag_a=protocol.define_liquid(name='Reagent A',display_color="#00E5FF")  # bright cyan
    bsa_reag_b=protocol.define_liquid(name='Reagent B',display_color="#FFD400")  # vivid yellow
    bsa_reag_c=protocol.define_liquid(name='Reagent C',display_color="#FF3D00")  # vivid orange-red
    biotin_azide=protocol.define_liquid(name='Biotin Azide',display_color="#FF4FBF")  # PINK (hot pink)
    copper_sulfate=protocol.define_liquid(name='CuSO4',display_color="#0066FF")  # BLUE (bright cobalt)
    tbta=protocol.define_liquid(name='TBTA',display_color="#7A00FF")  # vivid violet
    tcep_click=protocol.define_liquid(name='TCEP (Click)',display_color="#00C853")  # vivid green
    SP3_Beads=protocol.define_liquid(name='SP3 Beads',display_color="#111111")  # near-black (stands out)
    strep_beads=protocol.define_liquid(name='Streptavidin Beads',display_color="#795548")  # warm brown/gray (distinct from BSA)
    K2CO3=protocol.define_liquid(name='K2CO3',display_color="#FF6D00")  # bright amber
    tcep=protocol.define_liquid(name='TCEP',display_color="#00B8D4")  # teal (distinct from cyan above)
    empty_2mL=protocol.define_liquid(name='empty_2mL',display_color="#9E9E9E")  # gray
    IAA=protocol.define_liquid(name='IAA',display_color="#D50000")  # vivid red
    excess_lysis=protocol.define_liquid(name='Excess Lysis Buffer',display_color="#C51162")  # deep magenta
    epps_urea=protocol.define_liquid(name='2M Urea in EPPS',display_color="#00E676")  # neon green
    abs_ethanol=protocol.define_liquid(name='100% Ethanol',display_color="#1E88E5")  # blue (different shade from CuSO4)
    ethanol=protocol.define_liquid(name='Ethanol',display_color="#3949AB")  # indigo
    ethanol_80=protocol.define_liquid(name='80% Ethanol',display_color="#00ACC1")  # cyan-teal (unique)
    ACN=protocol.define_liquid(name='ACN',display_color="#E6FFF2")  # whitish green (light mint)
    epps=protocol.define_liquid(name='EPPS',display_color="#2979FF")  # bright azure (unique)
    trypsin=protocol.define_liquid(name='Trypsin in EPPS',display_color="#AEEA00")  # lime
    SDS=protocol.define_liquid(name='10% SDS',display_color="#263238")  # blue-gray / charcoal
    Urea6M_DPBS=protocol.define_liquid(name='6M Urea in DPBS',display_color="#E0F7FA")  # very light cyan
    DPBS=protocol.define_liquid(name='DPBS',display_color="#FFFFFF")  # white
    cacl2=protocol.define_liquid(name='CaCl2',display_color="#FF9100")  # orange (unique vs K2CO3)
    sample_liquids=[protocol.define_liquid(name=f'Sample {i+1}',display_color="#FFA000") for i in range(protocol.params.num_samples)]

    # Reservoir assignments for washes and digestion
    reservoir['A1'].load_liquid(liquid=bsa_reag_a, volume=15000)  
    reservoir['A2'].load_liquid(liquid=bsa_reag_b, volume=15000)  
    reservoir['A3'].load_liquid(liquid=bsa_reag_c, volume=15000)  
    reservoir['A4'].load_liquid(liquid=excess_lysis, volume=15000) 
    reservoir['A5'].load_liquid(liquid=SDS, volume = 15000)
    reservoir['A6'].load_liquid(liquid=Urea6M_DPBS, volume = 15000)
    reservoir['A7'].load_liquid(liquid=abs_ethanol,volume=15000)
    reservoir['A8'].load_liquid(liquid=ethanol_80, volume=15000)  # 80% Ethanol
    reservoir['A9'].load_liquid(liquid=DPBS, volume = 15000)
    reservoir['A10'].load_liquid(liquid=ACN, volume = 15000)
    reservoir['A11'].load_liquid(liquid=epps, volume=15000)  # 80% Ethanol for washes
    reservoir['A12'].load_liquid(liquid=epps_urea, volume=15000)  # 2M Urea in EPPS

    # Trypsin and CaCl2 for digestion
    temp_adapter['A1'].load_liquid(liquid=bsa_standard, volume=1000)  # 20 mg/ml BSA standard
    temp_adapter['A2'].load_liquid(liquid=copper_sulfate, volume=1500) #click
    temp_adapter['A3'].load_liquid(liquid=biotin_azide, volume=1500) #click
    temp_adapter['A4'].load_liquid(liquid=tcep_click, volume=1500) #click
    temp_adapter['A5'].load_liquid(liquid=tbta, volume=1500) #click
    temp_adapter['A6'].load_liquid(liquid=empty_2mL,volume=0)
    temp_adapter['D1'].load_liquid(liquid=SP3_Beads, volume=1000)  # Additional lysis buffer for SP3
    temp_adapter['D2'].load_liquid(liquid=tcep, volume=1000)  # 20 mg/ml BSA standard
    temp_adapter['D3'].load_liquid(liquid=K2CO3, volume=1000)  # Additional lysis buffer for SP3
    temp_adapter['D4'].load_liquid(liquid=IAA, volume=1000)  # 20 mg/ml BSA standard
    temp_adapter['D5'].load_liquid(liquid=cacl2, volume=500)  # CaCl2
    temp_adapter['D6'].load_liquid(liquid=trypsin, volume=2000)  # Trypsin in EPPS 
    temp_adapter['C6'].load_liquid(liquid=strep_beads, volume=2000)  # Trypsin in EPPS 

    # Load samples into the temp_adapter
    sample_rows = ['B', 'C', 'D', 'E']  # start at B
    temp_cols = 6  # adjust if your temp_adapter has fewer columns
    for i, liquid in enumerate(sample_liquids):
        row = sample_rows[i // temp_cols]
        col = (i % temp_cols) + 1
        well_name = f"{row}{col}"
        temp_adapter[well_name].load_liquid(liquid=liquid, volume=1500)  # adjust volume if needed

    # Load pipettes
    p50_multi = protocol.load_instrument('flex_8channel_50', 'left') 
    p1000_multi = protocol.load_instrument('flex_8channel_1000', 'right') 

    #Configure the p1000 pipette to use all channels
    p1000_multi.configure_nozzle_layout(style=ALL, tip_racks=[tips_200])

    # Steps 1: Add lysis buffer to column 1 of plate1. 
    p1000_multi.distribute(50, 
         reservoir['A4'],
         plate1[f'A{protocol.params.standards_col}'],
         rate = speed,
         mix_before=(1, 50),
         delay = 2,
         new_tip='once')

    #Step 3: Configure the p50 pipette to use single tip NOTE: this resets the pipettes tip racks! it doesn't
    p50_multi.configure_nozzle_layout(style=SINGLE, start="A1",tip_racks=[tips_50])

    # Step 2: move the 200uL tips to D4 and then the 50 uL partial tips to B3
    protocol.move_labware(labware=tips_200, new_location="B4", use_gripper=True)

    # Step 4: Transfer BSA standard (20 mg/ml) to first well of column 1
    p50_multi.transfer(50,
        temp_adapter['A1'],
        plate1[f'A{protocol.params.standards_col}'],
        rate = 0.35,
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
    # --- support up to 24 samples (B1–B6, C1–C6, D1–D6, E1–E6) ---
    if protocol.params.num_samples > 24:
        protocol.comment("num_samples > 24 requested; capping at 24.")
    num = min(protocol.params.num_samples, 24)

    # assign sample locations dynamically
    sample_locations = []
    for i in range(num):
        if i < 6:            # B1 to B6
            sample_locations.append(f'B{i + 1}')
        elif i < 12:         # C1 to C6
            sample_locations.append(f'C{i - 5}')
        elif i < 18:         # D1 to D6
            sample_locations.append(f'D{i - 11}')
        else:                # E1 to E6 (i = 18..23)
            sample_locations.append(f'E{i - 17}')

    # Create a list of rows that repeats based on protocol.params.num_samples
    row = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
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
                        disposal_vol=2)  # Distributing to three consecutive columns

    #Step 9: Load the p50 with full tip rack (don't need to)
    p50_multi.configure_nozzle_layout(style=ALL, tip_racks=[tips_50]) #, 

    #Step 10: Pipette triplicate of controls from plate1 column 1 to plate2 columns 1,2,3 
    p50_multi.distribute(5, 
                        plate1[f'A{protocol.params.standards_col}'], 
                        [plate2[f'A{i}'].bottom(z=0.1) for i in range(1, 4)],
                        rate= speed,
                        mix_before=(1, 10),
                        disposal_vol=5)

    #Step 12: Load the p1000 with full tip rack (don't need to)
    protocol.move_labware(labware=tips_1000, new_location='C3', use_gripper=True)

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
                        disposal_vol=5)

    #Step 16: move plate 2 to the heater shaker and incubate at 37c
    heater_shaker.open_labware_latch()
    protocol.move_labware(labware=plate2, new_location=heater_shaker, use_gripper=True)
    heater_shaker.set_and_wait_for_temperature(50)
    heater_shaker.close_labware_latch()
    heater_shaker.set_and_wait_for_shake_speed(500)
    protocol.delay(minutes=15)

    #Step 17 deactivate heater shaker and temp modules
    heater_shaker.deactivate_shaker()
    heater_shaker.deactivate_heater()
    heater_shaker.open_labware_latch()
    protocol.move_labware(labware=plate2, new_location='B2',use_gripper=True)

    # ---------------- Normalizing BCA Assay ----------------
    protocol.comment("Place BCA assay absorbance data in /var/lib/jupyter/notebooks/Data")

    # Pause the protocol until the user loads the file to /var/lib/jupyter/notebooks
    protocol.pause()

    # Tell the robot that new labware will be placed onto the deck
    protocol.move_labware(labware=plate1, new_location=protocol_api.OFF_DECK)
    protocol.move_labware(labware=plate2, new_location=protocol_api.OFF_DECK)
    protocol.move_labware(labware=plate3, new_location="B2", use_gripper=True)
    protocol.move_labware(labware=plate4, new_location=thermocycler)

    #Configure the p1000 pipette to use single tip NOTE: this resets the pipettes tip racks!
    protocol.move_labware(labware=tips_1000, new_location='C4', use_gripper=True)
    protocol.move_labware(labware=tips_200, new_location="B3", use_gripper=True)
    p1000_multi.configure_nozzle_layout(style=SINGLE, start="A1",tip_racks=[tips_200])

    # ---------------- Dummy normalization (10 samples) ----------------
    protocol.comment("Using dummy concentrations for normalization (skipping Excel + curve fit)")

    num_samples=10  # hardcoded test dataset

    # Dummy concentrations between 3–4 mg/mL
    dummy_conc_mg_ml=[3.2,3.4,3.6,3.8,3.3,3.5,3.7,3.9,3.25,3.75]

    final_vol_ul=protocol.params.final_volume
    target_conc_mg_ml=protocol.params.target_concentration

    # Calculate required sample volume to reach target concentration
    # sample_vol = (target * final_volume) / stock_concentration
    sample_vol_ul=[(target_conc_mg_ml*final_vol_ul)/c for c in dummy_conc_mg_ml]

    # Guardrails (cannot exceed final volume or be negative)
    sample_vol_ul=[min(max(v,0.0),final_vol_ul) for v in sample_vol_ul]

    # Calculate diluent volume
    diluent_vol_ul=[final_vol_ul - v for v in sample_vol_ul]

    # Build dataframe expected by downstream transfer loop
    normalized_samples=pd.DataFrame({
        "Sample":[f"Sample {i+1}" for i in range(num_samples)],
        "Protein Concentration (mg/mL)":dummy_conc_mg_ml,
        "Sample Volume (µL)":sample_vol_ul,
        "Diluent Volume (µL)":diluent_vol_ul
    })

    protocol.comment("Dummy normalization complete")

    rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    destination_wells  = [f'{rows[i % 8]}{(i // 8)+ 1}' for i in range(len(normalized_samples))]

    for i, row in normalized_samples.iterrows():
            source_well = sample_locations[i]
            normalized_volume = row['Sample Volume (µL)']
            diluent_volume = protocol.params.final_volume - normalized_volume
            destination_well = destination_wells[i]
            p1000_multi.transfer(diluent_volume, 
                                reservoir['A7'], 
                                plate3[destination_well].bottom(z=0.1), 
                                rate=0.5, 
                                new_tip='once')
            p1000_multi.transfer(normalized_volume, 
                                temp_adapter[source_well], 
                                plate3[destination_well].bottom(z=0.1), 
                                rate=0.5, 
                                new_tip='once')

    # ---------------- Click Reaction ----------------
    protocol.comment("Running click reaction")   

    #Pipette Biotin azide (A3), tbta (A5), cuso4 (A2), and tcep (A4)
    p1000_multi.transfer(5*(protocol.params.num_samples*2), 
                            temp_adapter['A3'], 
                            temp_adapter['A6'].bottom(z=0.1),
                            rate=speed,
                            mix_before=(1,10), 
                            disposal_vol=1,
                            new_tip='always')

    p1000_multi.transfer(30*(protocol.params.num_samples*2), 
                            temp_adapter['A5'], 
                            temp_adapter['A6'],
                            mix_before=(1,10),
                            rate=speed,
                            new_tip='always')

    p1000_multi.transfer(10*(protocol.params.num_samples*2), 
                            temp_adapter['A2'], 
                            temp_adapter['A6'], 
                            mix_before=(1,10),
                            new_tip='always')

    p1000_multi.transfer(10*(protocol.params.num_samples*2), 
                            temp_adapter['A4'], 
                            temp_adapter['A6'], 
                            new_tip='always')
    
    # Make sure the click reagents are well mixed
    click_volume = 6*(protocol.params.final_volume/50)

    # Track where tip racks are currently located
    tiprack_locations = {
        "tips_50": "A3",  # assume starts here
        "tips_200":"B3",    # assume starts her
        "tips_1000": "C4",   # assume starts here
    }

    def mix_add_click_reagents():
        volume_click_reaction = protocol.params.final_volume + click_volume
        location = temp_adapter['A6']
        pipette = None

        positions_mixing = [1, 1, 1]  # default fallback

        if volume_click_reaction < 50:
            positions_mixing = [1, 2, 3]
            if tiprack_locations["tips_50"] != "A3":
                protocol.move_labware(tips_50, new_location="A3", use_gripper=True)
                tiprack_locations["tips_50"] = "A3"
            if tiprack_locations["tips_200"] != "B4":
                protocol.move_labware(tips_200, new_location="B4", use_gripper=True)
                tiprack_locations["tips_200"] = "B4"
            p1000_multi.configure_nozzle_layout(style=SINGLE, start="A1", tip_racks=[tips_200])
            pipette = p1000_multi
            pipette = p50_multi

        elif 50 < volume_click_reaction < 200:
            positions_mixing = [1, 4, 9]
            if tiprack_locations["tips_1000"] != "C4":
                protocol.move_labware(tips_1000, new_location="C4", use_gripper=True)
                tiprack_locations["tips_1000"] = "C4"
            if tiprack_locations["tips_200"] != "B3":
                protocol.move_labware(tips_200, new_location="B3", use_gripper=True)
                tiprack_locations["tips_200"] = "B3"
            p1000_multi.configure_nozzle_layout(style=SINGLE, start="A1", tip_racks=[tips_200])
            pipette = p1000_multi

        elif 200 < volume_click_reaction < 1000:
            positions_mixing = [1, 6, 11]
            if tiprack_locations["tips_200"] != "B4":
                protocol.move_labware(tips_200, new_location="B4", use_gripper=True)
                tiprack_locations["tips_200"] = "B4"
            if tiprack_locations["tips_1000"] != "B3":
                protocol.move_labware(tips_1000, new_location="B3", use_gripper=True)
                tiprack_locations["tips_1000"] = "B3"
            p1000_multi.configure_nozzle_layout(style=SINGLE, start="A1", tip_racks=[tips_1000])
            pipette = p1000_multi

        else:
            pipette = p1000_multi

        # Perform mixing
        pipette.pick_up_tip()
        pipette.aspirate(protocol.params.final_volume / 2, location.bottom(z=positions_mixing[0]))
        pipette.dispense(protocol.params.final_volume / 2, location.bottom(z=positions_mixing[1]))
        pipette.aspirate(protocol.params.final_volume / 3, location.bottom(z=positions_mixing[2]))
        pipette.dispense(protocol.params.final_volume / 3, location.bottom(z=positions_mixing[0]))
        pipette.mix(3, protocol.params.final_volume, location.bottom(z=positions_mixing[0]))
        pipette.drop_tip()

        # Pipette the click reaction premix
        pipette.transfer(click_volume, 
                                temp_adapter['A6'], 
                                [plate3[i] for i in destination_wells],
                                rate=speed-0.1,
                                delay=2,
                                disposal_vol=0,
                                mix_before=(1, 6),
                                mix_after=(3,30),
                                new_tip='always')

    # Call the function
    mix_add_click_reagents()

    # Step 11: shake the sample plate for click reaction
    protocol.move_labware(labware=plate3, new_location=heater_shaker, use_gripper=True)
    heater_shaker.close_labware_latch()
    heater_shaker.set_and_wait_for_shake_speed(1000)
    protocol.delay(minutes=90)
    heater_shaker.deactivate_shaker()
    heater_shaker.open_labware_latch()
    protocol.move_labware(labware=plate3, new_location='B2', use_gripper=True)

    # ---------------- Protein Cleanup by SP3 ----------------
    protocol.comment("Desalting")
    
    # Move tip racks back to original locations in correct order
    if tiprack_locations["tips_50"] != "A3":
        protocol.move_labware(tips_50, new_location="A3", use_gripper=True)
        tiprack_locations["tips_50"] = "A3"
    if tiprack_locations["tips_200"] != "B4":
        protocol.move_labware(tips_200, new_location="B4", use_gripper=True)
        tiprack_locations["tips_200"] = "B4"
    if tiprack_locations["tips_1000"] != "B3":
        protocol.move_labware(tips_1000, new_location="B3", use_gripper=True)
        tiprack_locations["tips_1000"] = "B3"

    # mix the sp3 beads to homogenize. For a 10-plex sample, you need atleast 600 uL of beads.
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

    # Determine occupied columns on plate3 (block 0)
    cols_used, destination_columns, destination_wells_col, destination_wells_single = get_plate3_block(0)
    bead_volume = 50  # µL added per well in the SP3 bead step above
    starting_vol = protocol.params.final_volume + click_volume + bead_volume
    etoh_vol = (0.55 * starting_vol) / (1 - 0.55)  # µL of 100% EtOH to reach 55% v/v


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

        # After the final (3rd) wash, pull the last ~180 and residual with p50
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

    # resuspend in 6 M urea in DPBS with 0.2% SDS and move to B2
    protocol.move_labware(labware=plate3, new_location='B2', use_gripper=True)
    p1000_multi.transfer(500, 
                            reservoir['A6'], 
                            [well.bottom(z=10) for well in destination_wells_col],
                            mix_before=(1,150), 
                            mix_after=(3, 150), 
                            new_tip='once')

    #Step 9: Add 10% SDS with p50_multi to 0.2%  Note. 15 uL accounts for error when 
    p50_multi.transfer(15, 
                            reservoir['A5'], 
                            [well.bottom(z=1) for well in destination_wells_col],
                            rate=speed,
                            mix_before=(1,10), 
                            mix_after=(3, 10), 
                            new_tip='once')

    # Remove samples from beads and move to new wells on plate 3
    protocol.move_labware(labware=plate3, new_location=mag_block, use_gripper=True)
    protocol.delay(minutes=3)
    cols_used, shifted_destination_columns, shifted_destination_wells_col, shifted_destination_wells_single = get_plate3_block(1)

    for src, dst in zip(destination_wells_col, shifted_destination_wells_col):
        p1000_multi.transfer(
            510,
            src.bottom(z=1),
            dst.bottom(z=1),
            rate=0.25,
            new_tip='once'
        )
    protocol.move_labware(labware=plate3, new_location='B2', use_gripper=True)


    # ---------------- Reduction, Alkylation and Digestion ----------------
    #set the heater_shaker temp to 37 c for reduce
    heater_shaker.set_and_wait_for_temperature(37)

    #Configure the p1000 pipette to use single tip NOTE: this resets the pipettes tip racks!
    protocol.move_labware(labware=tips_1000_2, new_location="C4", use_gripper=True)
    p1000_multi.configure_nozzle_layout(style=SINGLE, start="A1",tip_racks=[tips_1000])

    #Pipette mix K2CO3, and TCEP and transfer to samples above the 200 uL liquid miniscus.
    Redu_mix = protocol.params.num_samples*50*1.2/2
    p1000_multi.transfer(Redu_mix, 
                        temp_adapter['A3'], 
                        temp_adapter['A6'], 
                        new_tip='once')

    p1000_multi.transfer(Redu_mix, 
                        temp_adapter['A4'], 
                        temp_adapter['A6'], 
                        new_tip='always')

    p1000_multi.distribute(50, 
                        temp_adapter['A6'], 
                        [w.top() for w in shifted_destination_wells_single], 
                        mix_before=(1,50),
                        disposal_vol=5,
                        new_tip='once')
    
    # Reduce for 30 minutes
    protocol.move_labware(labware=plate3, new_location=heater_shaker, use_gripper=True)
    heater_shaker.close_labware_latch()
    heater_shaker.set_and_wait_for_shake_speed(1000)
    protocol.delay(minutes=5) #################change me back to 30
    heater_shaker.deactivate_shaker()
    heater_shaker.open_labware_latch()

    # add IAA and alkylate for 30 minutes
    protocol.move_labware(labware=plate3, new_location='B2', use_gripper=True)
    p1000_multi.distribute(70, 
                            temp_adapter['A5'], 
                            [w.top() for w in shifted_destination_wells_single],
                            disposal_vol=5,
                            new_tip='once')

    protocol.move_labware(labware=plate3, new_location=heater_shaker, use_gripper=True)
    heater_shaker.close_labware_latch()
    heater_shaker.set_and_wait_for_shake_speed(1000)
    protocol.delay(minutes=5) #################change me back to 30
    heater_shaker.deactivate_shaker()
    heater_shaker.open_labware_latch()

    # move plate3 to B2 and add streptavidin magnetic beads then incubate for 2 hours on heater shaker
    protocol.move_labware(labware=plate3, new_location='B2', use_gripper=True)
    p1000_multi.distribute(100, 
                        temp_adapter['C6'], 
                        [w.bottom(z=15) for w in shifted_destination_wells_single],
                        mix_before=(3, 100), 
                        new_tip='always')

    protocol.move_labware(labware=plate3, new_location=heater_shaker, use_gripper=True)
    heater_shaker.close_labware_latch()
    heater_shaker.set_and_wait_for_shake_speed(1000)
    protocol.delay(minutes=90) #################change me back to 30
    heater_shaker.deactivate_shaker()
    heater_shaker.open_labware_latch()

    # ---------------- Washing the Streptavidin Beads ----------------
    #Move the full rack of p1000 tips
    protocol.move_labware(labware=tips_1000_2, new_location="C3", use_gripper=True)
    tiprack_locations["tips_1000_2"] = "C3"

    #Configure the p1000 pipette to use All tips 
    p1000_multi.configure_nozzle_layout(style=ALL, tip_racks=[tips_1000_2])
    protocol.move_labware(labware=plate3, new_location=mag_block, use_gripper=True)
    protocol.delay(minutes=15) ##this is a minimal time necessary to bind the beads in that much EtOH/lysis buffer.

    # Remove the supernatant from the beads leaving 200 uL in the bottom
    p1000_multi.pick_up_tip()
    for well in shifted_destination_wells_col:
        p1000_multi.aspirate(530, well.bottom(z=2), rate=0.1)
        protocol.delay(seconds=2)
        p1000_multi.blow_out(chute)
    p1000_multi.drop_tip()  

    # Remove remainder of EtOH and Add 0.2% SDS in DPBS and wash the beads three times
    p50_multi.configure_nozzle_layout(style=ALL, tip_racks=[tips_50])
    for i in range(3):
        # Ensure p1000 has a tip before any aspirate/distribute with new_tip='never'
        if not p1000_multi.has_tip:
            p1000_multi.pick_up_tip()

        # --- REMOVE SUPERNATANT (slow, above pellet) ---
        for well in shifted_destination_wells_col:
            p1000_multi.aspirate(180, well.bottom(z=1), rate=0.1)
            protocol.delay(seconds=2)
            p1000_multi.blow_out(chute)

        # --- Add 0.2% SDS in DPBS and wash the beads three times ---
        p1000_multi.distribute(
            180,
            reservoir['A9'], ################## Change me to correct washing buffers
            [w.top(z=-10) for w in shifted_destination_wells_col],
            rate=speed,
            new_tip='never'
        )
        p50_multi.distribute(4.4, 
                            reservoir['A5'], #
                            [well.top() for well in shifted_destination_wells_col],
                            disposal_vol=10,  
                            new_tip='once')

        # Shake, then back to magnet
        protocol.move_labware(labware=plate3, new_location=heater_shaker, use_gripper=True)
        heater_shaker.close_labware_latch()
        heater_shaker.set_and_wait_for_shake_speed(1000)
        protocol.delay(seconds=60)
        heater_shaker.deactivate_shaker()
        heater_shaker.open_labware_latch()
        protocol.move_labware(labware=plate3, new_location=mag_block, use_gripper=True)
        protocol.delay(minutes=3)

    # Remove remainder of SDS/DPBS and wash the beads three times with water
    p50_multi.configure_nozzle_layout(style=ALL, tip_racks=[tips_50])
    for i in range(3):
        # Ensure p1000 has a tip before any aspirate/distribute with new_tip='never'
        if not p1000_multi.has_tip:
            p1000_multi.pick_up_tip()

        # --- REMOVE SUPERNATANT (slow, above pellet) ---
        for well in shifted_destination_wells_col:
            p1000_multi.aspirate(180, well.bottom(z=1), rate=0.1)
            protocol.delay(seconds=2)
            p1000_multi.blow_out(chute)

        # --- Add 0.2% SDS in DPBS and wash the beads three times ---
        p1000_multi.distribute(
            180,
            reservoir['A9'], ################## Change me to correct washing buffers
            [w.top(z=-10) for w in shifted_destination_wells_col],
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

        # After the final (3rd) wash, pull the last ~180 and residual with p50
        if i == 2:
            if not p1000_multi.has_tip:
                p1000_multi.pick_up_tip()
            for well in shifted_destination_wells_col:
                p1000_multi.aspirate(180, well.bottom(z=1), rate=0.1)
                protocol.delay(seconds=2)
                p1000_multi.blow_out(chute)
            if p1000_multi.has_tip:
                p1000_multi.drop_tip()
            if not p50_multi.has_tip:
                p50_multi.pick_up_tip()
            for well in shifted_destination_wells_col:
                p50_multi.aspirate(50, well.bottom(z=0.1), rate=0.1)
                protocol.delay(seconds=2)
                p50_multi.blow_out(chute)
            p50_multi.drop_tip()
    
    # Resuspend in EPPS, Magnetize, and transfer to new column
    protocol.move_labware(labware=plate3, new_location=heater_shaker, use_gripper=True)
    heater_shaker.close_labware_latch()
    p1000_multi.distribute(500, 
                        reservoir['A11'], 
                        [well for well in shifted_destination_wells_col],
                        mix_before=(1,150), 
                        mix_after=(3, 150), 
                        new_tip='always')

    heater_shaker.set_and_wait_for_shake_speed(1000)
    protocol.delay(seconds=30)
    heater_shaker.deactivate_shaker()
    heater_shaker.open_labware_latch()
    protocol.move_labware(labware=plate3, new_location=mag_block, use_gripper=True)
    protocol.delay(minutes=3)

    # --- Next block (block 2) ---
    _, shift1_columns, shift1_wells_col, shift1_wells_single = get_plate3_block(2)

    # --- Safety check (96-well plate has 12 columns) ---
    if 3 * cols_used > len(plate3.columns()):
        raise ValueError(
            f"Not enough columns on plate3 for 3 blocks. "
            f"cols_used={cols_used}, need {3*cols_used} columns, plate has {len(plate3.columns())}."
    )
    for src, dst in zip(shifted_destination_wells_col, shift1_wells_col):
            p1000_multi.transfer(
                500,
                src.bottom(z=1),
                dst.bottom(z=1),
                rate=0.25,
                new_tip='once'
            )
    
    protocol.move_labware(labware=plate3, new_location='B2', use_gripper=True)
    #Configure the p1000 pipette to use single tip NOTE: this resets the pipettes tip racks!
    protocol.move_labware(labware=tips_1000_2, new_location='C4', use_gripper=True)
    p1000_multi.configure_nozzle_layout(style=SINGLE, start="A1",tip_racks=[tips_1000])

    # Add CaCl2, trypsin in epps, and move to shaker
    p1000_multi.distribute(2.5, 
                            temp_adapter['D5'], 
                            [plate3[w].top() for w in destination_wells],
                            disposal_vol=10, 
                            new_tip='once')

    p1000_multi.transfer(150, 
                            temp_adapter['D6'], 
                            [plate3[w].top() for w in destination_wells],
                            disposal_vol=0,
                            rate=speed,
                            mix_before=(1,150),
                            mix_after=(3,500),
                            new_tip='once')

    # Digest overnight
    heater_shaker.close_labware_latch()
    heater_shaker.set_and_wait_for_shake_speed(1000)
    protocol.delay(minutes=960)
    heater_shaker.deactivate_shaker()
    heater_shaker.open_labware_latch()
    protocol.comment("Samples have been digested")

    # Move digested samples to new columns again
    protocol.move_labware(labware=plate3, new_location=mag_block, use_gripper=True)
    protocol.delay(minutes=3)
    # --- Next block (block 3) ---
    _, shift2_columns, shift2_wells_col, shift2_wells_single = get_plate3_block(3)

    # --- Safety check (96-well plate has 12 columns) ---
    if 3 * cols_used > len(plate3.columns()):
        raise ValueError(
            f"Not enough columns on plate3 for 3 blocks. "
            f"cols_used={cols_used}, need {3*cols_used} columns, plate has {len(plate3.columns())}.")

    for src, dst in zip(shift1_wells_col, shift2_wells_col):
        p1000_multi.transfer(
            300,
            src.bottom(z=1),
            dst.bottom(z=1),
            rate=0.25,
            new_tip='once'
        )

    # ---------------- Post-Digestion TMT Reaction ----------------
    protocol.comment("Running TMT reaction")
    protocol.move_labware(labware=tmt_plate, new_location='A2')

    # plate3 needs to be accessible for ACN addition first
    protocol.move_labware(labware=plate3, new_location='B2', use_gripper=True)

    num_rows=8
    num_samples=protocol.params.num_samples
    cols_used=math.ceil(num_samples/num_rows)

    # --- plate3 block3 (shift2 block): columns [3*cols_used : 4*cols_used) ---
    block3_start_col_idx=3*cols_used
    block3_cols=plate3.columns()[block3_start_col_idx:block3_start_col_idx+cols_used]

    # Column anchors for multi-channel additions (A-row wells)
    block3_anchors=[col[0] for col in block3_cols]

    # Single-well list for SINGLE operations
    block3_sample_wells=[well for col in block3_cols for well in col][:num_samples]

    # --- Add ACN to 30% of 220 µL (66 µL total) using P1000 in style=ALL ---
    sample_vol=220
    acn_vol=0.30*sample_vol  # 66 µL
    acn_source=reservoir['A9']  # <-- set to your ACN well

    for dst_anchor in block3_anchors:
        p1000_multi.transfer(
            acn_vol,
            acn_source,
            dst_anchor.bottom(z=1.0),
            rate=0.5,
            new_tip='always',
            mix_after=(5, 150)
        )

    # move tips_200 (after ACN is added)
    protocol.move_labware(labware=tips_1000, new_location="C3", use_gripper=True)

    # --- Switch to P50 in SINGLE for TMT reagent addition ---
    p50_multi.configure_nozzle_layout(style=SINGLE, start="A1", tip_racks=[tips_50])

    tmt_sources=[tmt_plate.wells_by_name()[f"{protocol.params.tmt_row}{c}"] for c in range(1, num_samples+1)]
    tmt_add_vol=5  # or tmt_vol

    for source, dest in zip(tmt_sources, block3_sample_wells):
        p50_multi.transfer(
            tmt_add_vol,
            source,
            dest.bottom(z=0.5),
            disposal_vol=0,
            rate=0.2,
            new_tip='always',
            mix_after=(8, 40)
        )

    # React
    thermocycler.close_lid()
    protocol.delay(minutes=120)
    thermocycler.open_lid()

    # Quench with hydroxylamine
    for dest in block3_sample_wells:
        p50_multi.transfer(
            1,
            temp_adapter['D4'],
            dest.bottom(z=0.5),
            rate=speed,
            mix_before=(1, 5),
            mix_after=(2, 5),
            new_tip='always'
        )

    thermocycler.close_lid()
    protocol.delay(minutes=5)
    video_process.terminate()
    thermocycler.set_block_temperature(4)