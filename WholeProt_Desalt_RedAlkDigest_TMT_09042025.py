from opentrons import protocol_api
from opentrons.protocol_api import SINGLE, ALL
import datetime
import time

metadata = {
    'protocolName': 'TMT prep — start at EtOH removal on magnet',
    'author': 'Assistant',
    'description': 'Begins at the step removing EtOH from beads on magnet; includes everything needed to run from here.'
}

requirements = {
    "robotType": "Flex",
    "apiLevel": "2.21"
}

def add_parameters(parameters):
    parameters.add_int(
        variable_name="num_samples",
        display_name="Number of samples",
        description="Number of input samples",
        default=10,
        minimum=1,
        maximum=24
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
        description="Used to calculate the volume for 10 µg transfers",
        default=1.0,
        minimum=0.5,
        maximum=2.5,
        unit="µg/µL"
    )

def run(protocol: protocol_api.ProtocolContext):
    # -------- Minimal setup required to begin at the EtOH-removal-on-magnet step --------
    speed = 0.3

    # Modules
    heater_shaker = protocol.load_module('heaterShakerModuleV1', 'D1')
    thermocycler = protocol.load_module('thermocyclerModuleV2')
    temp_module = protocol.load_module('temperature module gen2', 'C1')
    mag_block = protocol.load_module('magneticBlockV1', 'D2')
    chute = protocol.load_waste_chute()

    # Labware
    plate3 = protocol.load_labware('nest_96_wellplate_2ml_deep', 'B2')  # sample plate
    plate4 = thermocycler.load_labware('opentrons_96_wellplate_200ul_pcr_full_skirt')  # labeling plate
    reservoir = protocol.load_labware('nest_12_reservoir_15ml', 'C2')

    # Tips
    partial_50 = protocol.load_labware('opentrons_flex_96_filtertiprack_50ul', 'A3')
    tips_1000 = protocol.load_labware('opentrons_flex_96_filtertiprack_1000ul', "B3")
    partial_1000 = protocol.load_labware('opentrons_flex_96_filtertiprack_1000ul', 'B4')

    # TMT source plate (bring on deck before labeling)
    tmt_plate = protocol.load_labware('opentrons_96_wellplate_200ul_pcr_full_skirt', location=protocol_api.OFF_DECK)

    # Temp adapter tubes (for hydroxylamine, CaCl2, trypsin)
    temp_adapter = temp_module.load_labware('opentrons_24_aluminumblock_nest_1.5ml_screwcap')
    heater_shaker.open_labware_latch()

    # Define liquids (names only for clarity; volumes optional here)
    hydroxylamine = protocol.define_liquid('hydroxylamine', '#8A2BE2')
    cacl2 = protocol.define_liquid('CaCl2', '#FF3300')
    trypsin = protocol.define_liquid('Trypsin in EPPS', '#9900FF')
    epps = protocol.define_liquid('200 mM EPPS', '#00FF99')
    epps_urea = protocol.define_liquid('2M Urea in EPPS', '#00FF77')

    # Optional: load minimal volumes for this stage
    temp_adapter['D4'].load_liquid(liquid=hydroxylamine, volume=200)
    temp_adapter['D5'].load_liquid(liquid=cacl2, volume=500)
    temp_adapter['D6'].load_liquid(liquid=trypsin, volume=2000)
    reservoir['A10'].load_liquid(liquid=protocol.define_liquid('80% EtOH', '#5A9BD5'), volume=15000)
    reservoir['A11'].load_liquid(liquid=epps, volume=15000)
    reservoir['A12'].load_liquid(liquid=epps_urea, volume=15000)

    # Instruments
    p50_multi = protocol.load_instrument('flex_8channel_50', 'left')
    p1000_multi = protocol.load_instrument('flex_8channel_1000', 'right')

    # Compute well targets
    rows = list("ABCDEFGH")
    destination_wells = [f'{rows[i % 8]}{(i // 8) + 1}' for i in range(protocol.params.num_samples)]
    num_full_columns = (protocol.params.num_samples + 7) // 8
    destination_columns = plate3.columns()[:num_full_columns]
    destination_wells_col = [col[0] for col in destination_columns]  # A1, A2, ...

    # Ensure plate is on magnet for this starting step
    protocol.move_labware(labware=plate3, new_location=mag_block, use_gripper=True)
    protocol.delay(minutes=3)

    # ----------------------- YOUR REQUESTED START POINT -----------------------
    # Remove the EtOH from the beads leaving 200 uL in the bottom
    p1000_multi.configure_nozzle_layout(style=ALL, tip_racks=[tips_1000])
    p1000_multi.pick_up_tip()
    for well in destination_wells_col:
        p1000_multi.aspirate(900, well.bottom(z=2), rate=0.25)
        p1000_multi.blow_out(chute)
    # -------------------------------------------------------------------------

    # Remove remainder of EtOH and wash beads 3× with 80% EtOH (200 µL)
    for i in range(3):
        for well in destination_wells_col:
            p1000_multi.aspirate(200, well.bottom(z=2), rate=speed-0.1)
            p1000_multi.blow_out(chute)

        # Off magnet, add 80% EtOH, brief shake, then back to magnet
        protocol.move_labware(labware=plate3, new_location=heater_shaker, use_gripper=True)
        heater_shaker.close_labware_latch()
        p1000_multi.distribute(
            200,
            reservoir['A10'],
            [w.top() for w in destination_wells_col],
            rate=0.35,
            new_tip='never'
        )
        heater_shaker.set_and_wait_for_shake_speed(1000)
        protocol.delay(seconds=60)
        heater_shaker.deactivate_shaker()
        heater_shaker.open_labware_latch()

        protocol.move_labware(labware=plate3, new_location=mag_block, use_gripper=True)
        protocol.delay(minutes=3)

        # After the final (3rd) wash, remove the last ~200 µL
        if i == 2:
            for well in destination_wells_col:
                p1000_multi.aspirate(200, well.bottom(z=2), rate=speed-0.1)
                p1000_multi.blow_out(chute)
    p1000_multi.drop_tip()
    # Resuspend in 2 M urea in EPPS (147.5 µL) for digestion setup
    protocol.move_labware(labware=plate3, new_location='B2', use_gripper=True)
    p1000_multi.transfer(
        147.5,
        reservoir['A12'],
        [well for well in destination_wells_col],
        mix_before=(1, 150),
        mix_after=(3, 150),
        new_tip='always'
    )

    # Switch to SINGLE (fresh 1000 rack if desired; using same tips_1000 here)
    p1000_multi.configure_nozzle_layout(style=SINGLE, start="A1", tip_racks=[tips_1000])

    # Add CaCl2 (2.5 µL) then trypsin (150 µL)
    p1000_multi.distribute(
        2.5,
        temp_adapter['D5'],
        [plate3[w].top() for w in destination_wells],
        disposal_vol=10,
        new_tip='once'
    )
    p1000_multi.transfer(
        150,
        temp_adapter['D6'],
        [plate3[w].top() for w in destination_wells],
        disposal_vol=0,
        rate=speed,
        mix_before=(1, 150),
        mix_after=(3, 500),
        new_tip='once'
    )

    # Overnight digestion on shaker
    protocol.move_labware(labware=plate3, new_location=heater_shaker, use_gripper=True)
    heater_shaker.close_labware_latch()
    heater_shaker.set_and_wait_for_shake_speed(1000)
    protocol.delay(minutes=960)
    heater_shaker.deactivate_shaker()
    heater_shaker.open_labware_latch()
    protocol.comment("Digestion complete.")

    