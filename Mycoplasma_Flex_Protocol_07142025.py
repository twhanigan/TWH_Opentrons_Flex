from opentrons import protocol_api
from opentrons.protocol_api import SINGLE, ALL
import subprocess

metadata = {
    'protocolName': 'Mycoplasma Detection PCR Protocol Tube-based',
    'author': 'Assistant',
    'description': 'Automated PCR setup and gel preparation for Mycoplasma detection.',
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
        variable_name="volume_sample",
        display_name="volume sample",
        description="How much to aspirate from each sample.",
        default=2.5,
        minimum=1,
        maximum=5,
        unit="µL"
    )
def run(protocol: protocol_api.ProtocolContext):

    # variables for run
    numtotalSamples = protocol.params.num_samples*1.5 + (2*protocol.params.num_replicates)
    mastermix_vol = 12.5
    primer_vol = 1
    water_vol = 9
    # Step 4: Gel preparation and loading (manual step for now)
    protocol.comment("This protcol runs pcr assay for mycoplasma contamination.")
    
    #Start recording the video
    video_process = subprocess.Popen(["python3", "/var/lib/jupyter/notebooks/record_video_myco.py"])
    
    # Load modules
    heater_shaker = protocol.load_module('heaterShakerModuleV1', 'D1')
    thermocycler = protocol.load_module('thermocyclerModuleV2')
    temp_module = protocol.load_module('temperature module gen2', 'C1')
    mag_block = protocol.load_module('magneticBlockV1', 'D2')
    chute = protocol.load_waste_chute()
    temp_adapter = temp_module.load_labware('opentrons_24_aluminumblock_nest_1.5ml_screwcap')
    reservoir = protocol.load_labware('nest_12_reservoir_15ml', 'C2')

    # Load labware
    pcr_plate = thermocycler.load_labware('nest_96_wellplate_100ul_pcr_full_skirt')
    partial_50 = protocol.load_labware(load_name="opentrons_flex_96_tiprack_50ul",location="A2")
    tips_50 = protocol.load_labware('opentrons_flex_96_tiprack_50ul', location='A3')
    tips_1000 = protocol.load_labware(load_name="opentrons_flex_96_filtertiprack_1000ul",location="B3")
    
    # Load the liquids
    mastermix = protocol.define_liquid(name='mastermix', display_color="#FF00FF")  # Brown
    primer_mix = protocol.define_liquid(name = 'primer_mix', display_color="#00FFFF")
    nuclease_free_water = protocol.define_liquid(name = 'nuclease_free_water', display_color="#FFFF00")
    positive_control = protocol.define_liquid(name = 'positive_control', display_color="#FF3300")
    neg_control = protocol.define_liquid(name = 'neg_control', display_color="#FF5E00")
    empty_epp = protocol.define_liquid(name = 'empty_eppendorf', display_color="#000011")
    loading_buff = protocol.define_liquid(name = 'loading_buff', display_color="#220011")
    sample_tubes = [protocol.define_liquid(name = f'Sample {i + 1}', display_color="#FFA000",) for i in range(protocol.params.num_samples)]


    # Assign sample/liquid positions
    temp_adapter['A1'].load_liquid(liquid=mastermix, volume=1000)  # 20 mg/ml BSA standard
    temp_adapter['A2'].load_liquid(liquid=primer_mix, volume=1000)  # Additional lysis buffer for SP3
    temp_adapter['A3'].load_liquid(liquid=nuclease_free_water, volume=1000)  # 20 mg/ml BSA standard
    temp_adapter['A4'].load_liquid(liquid=positive_control, volume=1000)  # Additional lysis buffer for SP3
    temp_adapter['A5'].load_liquid(liquid=neg_control, volume=1000)  # 20 mg/ml BSA standard
    temp_adapter['A6'].load_liquid(liquid=empty_epp, volume=1000)  # 20 mg/ml BSA standard
    reservoir['A2'].load_liquid(liquid=loading_buff, volume=5000)
    
    #open the thermocycler lid
    thermocycler.open_lid()
    temp_module.set_temperature(celsius=10)
    
    # Load pipettes and configure nozzle layout
    p50_multi = protocol.load_instrument('flex_8channel_50', 'left') #, tip_racks=[tips_50]
    p1000_multi = protocol.load_instrument('flex_8channel_1000', 'right',tip_racks=[tips_1000]) #, 
    p50_multi.configure_nozzle_layout(style=SINGLE, start="A1",tip_racks=[partial_50])
    p1000_multi.configure_nozzle_layout(style=SINGLE, start="A1",tip_racks=[tips_1000])

    # assign sample locations dynamically
    sample_locations = []
    row = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    for i in range(protocol.params.num_samples):
        if i < 6:  # B1 to B6
            sample_locations.append(f'B{i + 1}')
        elif i < 12:  # C1 to C6
            sample_locations.append(f'C{i - 5}')
        elif i < 18:  # D1 to D6
            sample_locations.append(f'D{i - 11}')
        elif i < 23:
            sample_locations.append(f'A{i - 16}')  # Stop if we exceed the number of available rows/columns
        else:
            print('Too many samples')
    
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
    
    # cool the thermocycler
    thermocycler.set_block_temperature(4)

    #Add the positive/negative control
    p50_multi.distribute(protocol.params.volume_sample,
                         temp_adapter['A4'],
                         [pcr_plate[well].bottom(z=0.1) for well in sample_well_map["positive_control"]],
                         rate=0.5,
                         mix_before=(1, 10))

    # Transfer negative control to B1–B3
    p50_multi.distribute(protocol.params.volume_sample,
                         temp_adapter['A5'],
                         [pcr_plate[well].bottom(z=0.1) for well in sample_well_map["neg_control"]],
                         rate=0.5,
                         mix_before=(1, 10))
    
    # Transfer samples
    for sample_idx in range(protocol.params.num_samples):
        tube_loc = sample_locations[sample_idx]
        wells = sample_well_map[f"sample_{sample_idx+1}"]
        p50_multi.distribute(
            protocol.params.volume_sample,
            temp_adapter[tube_loc],
            [pcr_plate[well].bottom(z=0.1) for well in wells],
            rate=0.5,
            mix_before=(1, 5),
            disposal_vol=1)

    # Step 1: Distribute mastermix, primer mix, and water into PCR plate starting at row B1:
    p1000_multi.distribute((mastermix_vol*numtotalSamples*protocol.params.num_replicates), 
                            temp_adapter['A1'].bottom(z=2), 
                            temp_adapter['A6'], 
                            new_tip='once',
                            mix_before=(1, 10),
                            rate= 0.3)

    p1000_multi.distribute((primer_vol*numtotalSamples*protocol.params.num_replicates), 
                            temp_adapter['A2'].bottom(z=2), 
                            temp_adapter['A6'], 
                            new_tip='once')

    p1000_multi.distribute((water_vol*numtotalSamples*protocol.params.num_replicates), 
                            temp_adapter['A3'].bottom(z=2), 
                            temp_adapter['A6'], 
                            mix_after=(3, 40),
                            new_tip='once')

    mastermix_wells = [well for key, wells in sample_well_map.items() for well in wells]
    # distribute mastermix
    p50_multi.distribute(
        22.5,
        temp_adapter['A6'],
        [pcr_plate[well].bottom(z=5) for well in mastermix_wells],
        new_tip='once',
        disposal_vol=5,
        rate=0.5
    )

    # Step 3: Run thermocycling conditions
    thermocycler.set_lid_temperature(105)
    thermocycler.close_lid()

    # Stop video recording after the main task is completed
    video_process.terminate()
    
    protocol.comment('Running thermocycler for 30 minutes')
    # Initial Denaturation
    thermocycler.execute_profile(
        steps=[
            {'temperature': 94, 'hold_time_seconds': 120},
        ],
        repetitions=1,
        block_max_volume=25
    )

    # PCR Cycles (34x)
    thermocycler.execute_profile(
        steps=[
            {'temperature': 94, 'hold_time_seconds': 30},   # Denaturation
            {'temperature': 55, 'hold_time_seconds': 120},  # Annealing
            {'temperature': 72, 'hold_time_seconds': 60},   # Extension
        ],
        repetitions=34,
        block_max_volume=25
    )

    # Final Extension
    thermocycler.execute_profile(
        steps=[
            {'temperature': 72, 'hold_time_seconds': 60},
        ],
        repetitions=1,
        block_max_volume=25
    )

    thermocycler.deactivate_lid()
    thermocycler.set_block_temperature(4)  # Hold at 4°C
    #thermocycler.open_lid()

    # Step 4: Gel preparation and loading (manual step for now)
    protocol.comment("After PCR, analyze products on a 2% agarose gel stained with ethidium bromide")