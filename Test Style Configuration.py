from opentrons import protocol_api
from opentrons.protocol_api import SINGLE, ALL
from opentrons import types

metadata = {
    'protocolName': 'Edited Practice Next Tip',
    'author': 'Assistant',
    'description': 'Serial dilution of BSA standard and sample processing. This includes cooling samples to 4c, heating plate to 37c with shaking and recording a video of the whole process. Place BSA Standard in A1, Lysis buffer in A2, change the number of samples and place samples in row B starting at B1. MINIMUM Sample volumen in eppendorf tubes is 40 uL. '
}

requirements = {
    "robotType": "Flex",
    "apiLevel": "2.21"
}

def run(protocol: protocol_api.ProtocolContext):
    chute = protocol.load_waste_chute()
    partial_50 = protocol.load_labware(load_name="opentrons_flex_96_filtertiprack_50ul",location="B3")

    # Load pipettes
    p50_multi = protocol.load_instrument('flex_8channel_50', 'left')

    #Test
    p50_multi.configure_nozzle_layout(style=ALL, tip_racks=[partial_50])
    p50_multi.pick_up_tip()
    p50_multi.drop_tip()

    p50_multi.configure_nozzle_layout(style=SINGLE, start="A1", tip_racks=[partial_50])
    p50_multi.pick_up_tip()
    p50_multi.drop_tip()

    p50_multi.configure_nozzle_layout(style=ALL, tip_racks=[partial_50])
    p50_multi.pick_up_tip()
    p50_multi.drop_tip()
