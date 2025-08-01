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

def run(protocol: protocol_api.ProtocolContext):

    # Enter the number of samples 
    speed= 0.5
    target_concentration = 1
    num_samples = 1
    num_replicates = 2
    numtotalSamples = num_samples + (2*num_replicates)
    reaction_vol = 25
    mastermix_vol = 12.5
    primer_vol = 1
    sample_vol = 2.5
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
    empty_epp = protocol.define_liquid(name = 'neg_control', display_color="#000011")
    loading_buff = protocol.define_liquid(name = 'loading_buff', display_color="#220011")
    sample_tubes = protocol.define_liquid(name=samples, display_color= "#230442")

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
    
    # Load pipettes
    p50_multi = protocol.load_instrument('flex_8channel_50', 'left') #, tip_racks=[tips_50]
    p1000_multi = protocol.load_instrument('flex_8channel_1000', 'right',tip_racks=[tips_1000]) #, 

    #Configure the p50 pipette to use single tip NOTE: this resets the pipettes tip racks!
    p50_multi.configure_nozzle_layout(style=SINGLE, start="A1",tip_racks=[partial_50])
    p1000_multi.configure_nozzle_layout(style=SINGLE, start="A1",tip_racks=[tips_1000])

    # assign sample locations dynamically
    sample_locations = []
    row = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    for i in range(num_samples):
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
    def get_next_wells(start_index, num_replicates, used_wells):
        col = 1 + (start_index // (num_rows // num_replicates))
        row_pair_index = start_index % (num_rows // num_replicates)
        wells = []
        for i in range(num_replicates):
            row_letter = row_letters[(row_pair_index * num_replicates) + i]
            well = f"{row_letter}{col}"
            wells.append(well)
            used_wells.add(well)
        return wells

    # create mapping of samples and controls
    used_wells = set()
    sample_well_map = {}

    # allocate wells for positive control
    sample_well_map["positive_control"] = get_next_wells(0, num_replicates, used_wells)

    # allocate wells for negative control
    sample_well_map["neg_control"] = get_next_wells(1, num_replicates, used_wells)

    # allocate wells for each sample
    for sample_idx in range(num_samples):
        next_wells = get_next_wells(sample_idx + 2, num_replicates, used_wells)
        sample_well_map[f"sample_{sample_idx+1}"] = next_wells
    
    #Add the positive control and no template control to the number of samples
    # Transfer positive control to A1
    p50_multi.distribute(2.5,
                         temp_adapter['A4'],
                         [pcr_plate[well].bottom(z=0.1) for well in sample_well_map["positive_control"]],
                         rate=speed,
                         mix_before=(1, 10),
                         disposal_vol=5)

    # Transfer negative control to B1–B3
    p50_multi.distribute(2.5,
                         temp_adapter['A5'],
                         [pcr_plate[well].bottom(z=0.1) for well in sample_well_map["neg_control"]],
                         rate=speed,
                         mix_before=(1, 10),
                         disposal_vol=5)
    
    # Transfer samples
    for sample_idx in range(num_samples):
        tube_loc = sample_locations[sample_idx]
        wells = sample_well_map[f"sample_{sample_idx+1}"]

        p50_multi.distribute(
            2.5,
            temp_adapter[tube_loc],
            [pcr_plate[well].bottom(z=0.1) for well in wells],
            rate=speed,
            mix_before=(1, 10),
            disposal_vol=5)

    # Step 1: Distribute mastermix, primer mix, and water into PCR plate starting at row B1:
    p1000_multi.distribute((mastermix_vol*numtotalSamples*num_replicates), 
                            temp_adapter['A1'].bottom(z=2), 
                            temp_adapter['A6'], 
                            new_tip='once',
                            mix_before=(1, 10),
                            rate= speed-0.2)

    p1000_multi.distribute((primer_vol*numtotalSamples*num_replicates), 
                            temp_adapter['A2'].bottom(z=2), 
                            temp_adapter['A6'], 
                            new_tip='once')

    p1000_multi.distribute((water_vol*numtotalSamples*num_replicates), 
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
        rate=speed
    )

    # Step 3: Run thermocycling conditions
    thermocycler.close_lid()
    thermocycler.set_lid_temperature(105)
    protocol.comment('Running thermocycler for 30 minutes')
    thermocycler.set_block_temperature(4, hold_time_minutes=1)
    thermocycler.execute_profile(
        steps=[
            {'temperature': 95, 'hold_time_seconds': 180},  # Initial Denaturation
        ] + [
            {'temperature': 95, 'hold_time_seconds': 30},   # Denaturation
            {'temperature': 55, 'hold_time_seconds': 30},   # Annealing
            {'temperature': 72, 'hold_time_seconds': 30},   # Extension
        ] * 35
        + [
            {'temperature': 72, 'hold_time_seconds': 300}   # Final Extension, 5 min
        ],
        repetitions=34
    )

    thermocycler.set_block_temperature(4)  # Hold at 4°C
    #thermocycler.open_lid()

    # Stop video recording after the main task is completed
    video_process.terminate()
    
    # Step 4: Gel preparation and loading (manual step for now)
    protocol.comment("After PCR, analyze products on a 2% agarose gel stained with ethidium bromide")
