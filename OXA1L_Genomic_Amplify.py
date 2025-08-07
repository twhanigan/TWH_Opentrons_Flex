from opentrons import protocol_api
from opentrons.protocol_api import SINGLE, ALL
import subprocess

metadata = {
    'protocolName': 'OXA1L PCR Reaction',
    'author': 'Assistant',
    'description': 'PCR Amplification of OXA1L'
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
        variable_name="num_replicates",
        display_name="Number of replicates",
        description="Number of replicates for each sample",
        default=1,
        choices=[
            {"display_name": "1", "value": 1},
            {"display_name": "2", "value": 2},
            {"display_name": "3", "value": 3},
        ]
    )
    parameters.add_float(
        variable_name="reaction_vol",
        display_name="reaction volume",
        description="Volume of phusion reaction",
        default=50,
        minimum=1,
        maximum=200,
        unit="µL"
    )
    parameters.add_float(
        variable_name="sample_conc",
        display_name="sample concentration",
        description="Concentration of samples",
        default=100,
        minimum=1,
        maximum=1000,
        unit="ng/µL"
    )

def run(protocol: protocol_api.ProtocolContext):
    protocol.comment(
        f"Running OXA1L PCR Reaction on {protocol.params.num_samples}")
    
    # Change these if not using 96-well
    numtotalSamples = protocol.params.num_samples + 2
    sample_vol  = 250/protocol.params.sample_conc
    buffer_vol = 10
    primer_vol = 2.5
    dmso_vol = 1
    phusion_vol = 0.5
    speed= 0.35
    dntp_vol = 1
    mm_vol = protocol.params.reaction_vol -sample_vol
    water_vol = protocol.params.reaction_vol -(sample_vol+buffer_vol+dntp_vol+(primer_vol*2)+dmso_vol+phusion_vol)

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
    thermocycler.open_lid()
    thermocycler.set_block_temperature(4)  # Hold at 4°C
    
    # Load labware
    tips_50 = protocol.load_labware('opentrons_flex_96_filtertiprack_50ul', 'A2')
    tips_200 = protocol.load_labware('opentrons_flex_96_filtertiprack_200ul', 'A3')
    tips_1000 = protocol.load_labware(load_name="opentrons_flex_96_filtertiprack_1000ul",location="B3")
    reservoir = protocol.load_labware('nest_12_reservoir_15ml', 'C2')

    # Load PCR plate on thermocycler
    pcr_plate = thermocycler.load_labware('opentrons_96_wellplate_200ul_pcr_full_skirt')
    
    # Liquid definitionss
    OXA1L_F_Primer = protocol.define_liquid(name = 'OXA1L_F_Primer', display_color="#704848",)
    OXA1L_R_Primer = protocol.define_liquid(name = 'OXA1L_R_Primer', display_color="#704300",)
    HF_Buffer = protocol.define_liquid(name = 'HF_Buffer', display_color="#704900",)
    Phusion = protocol.define_liquid(name = 'Phusion', display_color="#701100",)
    ddH2O = protocol.define_liquid(name = 'ddH2O', display_color="#704900",)
    DMSO = protocol.define_liquid(name = 'DMSO', display_color="#704300",)
    dNTPs = protocol.define_liquid(name='dNTPs', display_color="#FFC0CB")  # Pink
    positive_control = protocol.define_liquid(name = 'positive_control', display_color="#FF3300")
    neg_control = protocol.define_liquid(name = 'neg_control', display_color="#FF5E00")
    empty_epp = protocol.define_liquid(name = 'empty_eppendorf', display_color="#000011")
    sample_tubes = [protocol.define_liquid(name = f'Sample {i + 1}', display_color="#FFA000",) for i in range(protocol.params.num_samples)]
    
    temp_adapter['A1'].load_liquid(liquid=ddH2O, volume=1500) #click
    temp_adapter['A2'].load_liquid(liquid=HF_Buffer, volume=1500) #click
    temp_adapter['A3'].load_liquid(liquid=dNTPs, volume=1500) #click
    temp_adapter['A4'].load_liquid(liquid=OXA1L_F_Primer, volume=1500) #click
    temp_adapter['A5'].load_liquid(liquid=OXA1L_R_Primer, volume=1500) #click
    temp_adapter['A6'].load_liquid(liquid=DMSO, volume=1500) #click
    temp_adapter['B1'].load_liquid(liquid=Phusion, volume=1500) #click
    temp_adapter['B2'].load_liquid(liquid=positive_control, volume=1000)  # Additional lysis buffer for SP3
    temp_adapter['B3'].load_liquid(liquid=neg_control, volume=1000) 
    temp_adapter['B4'].load_liquid(liquid=empty_epp, volume=1000) 

    # Load pipettes
    p50_multi = protocol.load_instrument('flex_8channel_50', 'left')
    p1000_multi = protocol.load_instrument('flex_8channel_1000', 'right')
    p50_multi.configure_nozzle_layout(style=SINGLE, start="A1",tip_racks=[tips_50])
    p1000_multi.configure_nozzle_layout(style=SINGLE, start="A1",tip_racks=[tips_1000])
    
    # assign sample locations dynamically
    if protocol.params.num_samples > 12:
        raise ValueError("Too many samples: only 12 wells available in rows C and D (C1–C6 and D1–D6).")
    sample_locations = []
    for i in range(protocol.params.num_samples):
            row = 'C' if i < 6 else 'D'
            col = (i % 6) + 1
            sample_locations.append(f"{row}{col}")    

    # define row letters in order
    samples = [temp_adapter[i] for i in sample_locations]
    row_letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    num_rows = len(row_letters)
    num_columns = 12

    # generate well positions for replicates placed vertically
    # start from column 1, row pairs A/B, C/D, etc.
    def get_next_wells(start_index, number_replicates, used_wells):
        col = 1 + (start_index // (num_rows // number_replicates))
        row_pair_index = start_index % (num_rows // number_replicates)
        wells = []
        for i in range(number_replicates):
            row_letter = row_letters[(row_pair_index * number_replicates) + i]
            well = f"{row_letter}{col}"
            wells.append(well)
            used_wells.add(well)
        return wells

    # create mapping of samples and controls
    used_wells = set()
    sample_well_map = {}

    # allocate wells for positive control
    sample_well_map["positive_control"] = get_next_wells(0, protocol.params.num_replicates, used_wells)

    # allocate wells for negative control
    sample_well_map["neg_control"] = get_next_wells(1, protocol.params.num_replicates, used_wells)

    # allocate wells for each sample
    for sample_idx in range(protocol.params.num_samples):
        next_wells = get_next_wells(sample_idx + 2, protocol.params.num_replicates, used_wells)
        sample_well_map[f"sample_{sample_idx+1}"] = next_wells
    
    thermocycler.set_block_temperature(4)
    #Add the positive control and no template control to the number of samples
    # Transfer positive control to A1
    p50_multi.distribute(sample_vol,
                         temp_adapter['B2'],
                         [pcr_plate[well].bottom(z=0.1) for well in sample_well_map["positive_control"]],
                         rate=0.5,
                         mix_before=(1, 10))

    # Transfer negative control to B1–B3
    p50_multi.distribute(sample_vol,
                         temp_adapter['B3'],
                         [pcr_plate[well].bottom(z=0.1) for well in sample_well_map["neg_control"]],
                         rate=0.5,
                         mix_before=(1, 10))
    
    # Transfer samples
    for sample_idx in range(protocol.params.num_samples):
        tube_loc = sample_locations[sample_idx]
        wells = sample_well_map[f"sample_{sample_idx+1}"]
        p50_multi.distribute(
            sample_vol,
            temp_adapter[tube_loc],
            [pcr_plate[well].bottom(z=0.1) for well in wells],
            rate=0.5,
            mix_before=(1, 5),
            disposal_vol=1)
    
    # Make mastermix in empty eppendorf:
    p1000_multi.distribute((water_vol*numtotalSamples*protocol.params.num_replicates), 
                            temp_adapter['A1'], 
                            temp_adapter['B4'],
                            rate=0.5, 
                            new_tip='always')

    p1000_multi.distribute((buffer_vol*numtotalSamples*protocol.params.num_replicates), 
                            temp_adapter['A2'], 
                            temp_adapter['B4'], 
                            new_tip='always',
                            mix_before=(1, 20),
                            rate=speed-0.1,
                            delay=2)
   
    p50_multi.distribute((dntp_vol*numtotalSamples*protocol.params.num_replicates), 
                            temp_adapter['A3'], 
                            temp_adapter['B4'], 
                            new_tip='always',
                            mix_before=(1, 10),
                            rate=0.5)

    p50_multi.distribute((primer_vol*numtotalSamples*protocol.params.num_replicates), 
                            temp_adapter['A4'], 
                            temp_adapter['B4'], 
                            rate=0.5,
                            mix_before=(1,10),
                            new_tip='always')

    p50_multi.distribute((primer_vol*numtotalSamples*protocol.params.num_replicates), 
                            temp_adapter['A5'], 
                            temp_adapter['B4'],
                            rate=0.5,
                            mix_before=(1,10), 
                            new_tip='always')
    
    p50_multi.distribute((dmso_vol*numtotalSamples*protocol.params.num_replicates), 
                            temp_adapter['A6'], 
                            temp_adapter['B4'],
                            mix_before=(1, 10), 
                            new_tip='always',
                            rate=speed)

    p50_multi.distribute((phusion_vol*numtotalSamples*protocol.params.num_replicates), 
                            temp_adapter['B1'], 
                            temp_adapter['B4'], 
                            mix_before=(1,3),
                            rate = speed-0.1,
                            delay = 2,
                            new_tip='always')


    # Transfer mastermix
    mastermix_wells = [well for key, wells in sample_well_map.items() for well in wells]
    p1000_multi.distribute(
        mm_vol,
        temp_adapter['B4'],
        [pcr_plate[well].bottom(z=5) for well in mastermix_wells],
        new_tip='always',
        disposal_vol=5,
        rate=0.5
    )
    # Close thermocycler lid and set temperature
    thermocycler.set_lid_temperature(105)
    thermocycler.close_lid()

    # 7. Run PCR cycles
    # Initial denaturation
    thermocycler.execute_profile(
        steps=[{'temperature': 98, 'hold_time_seconds': 30}],
        repetitions=1,
        block_max_volume=50
    )
    
    # 34 cycles of PCR
    thermocycler.execute_profile(
        steps=[
            {'temperature': 98, 'hold_time_seconds': 10},
            {'temperature': 60, 'hold_time_seconds': 15},
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
    
    # Stop video recording after the main task is completed
    video_process.terminate()
    
    # Hold at 4°C
    thermocycler.set_block_temperature(4)
    
    # Step 4: Gel preparation and loading (manual step for now)
    protocol.comment("After PCR, analyze products on a 0.75% agarose gel stained with ethidium bromide")