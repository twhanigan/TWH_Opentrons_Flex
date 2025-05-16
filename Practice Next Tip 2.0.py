from opentrons import protocol_api
from opentrons.protocol_api import SINGLE, PARTIAL_COLUMN, ALL
from opentrons import types

metadata = {
    'protocolName': 'Practice Next Tip',
    'author': 'Assistant',
    'description': 'Serial dilution of BSA standard and sample processing. This includes cooling samples to 4c, heating plate to 37c with shaking and recording a video of the whole process. Place BSA Standard in A1, Lysis buffer in A2, change the number of samples and place samples in row B starting at B1. MINIMUM Sample volumen in eppendorf tubes is 40 uL. '
}

requirements = {
    "robotType": "Flex",
    "apiLevel": "2.21"
}

def run(protocol: protocol_api.ProtocolContext):
    chute = protocol.load_waste_chute()
    partial_50 = protocol.load_labware(load_name="opentrons_flex_96_filtertiprack_50ul",location="A2")
    
    # Load pipettes
    p50_multi = protocol.load_instrument('flex_8channel_50', 'left')

    def search_tip(channels):
        if remaining_tips != []:
            for i in remaining_tips:
                if len(i) >= channels:
                    tip = i
                    remaining_tips.pop(i)
                    break
            return tip
        else:
            raise Exception("No tips available")

    def pick(channels, pip = p50_multi):
        if 2 <= channels <=7:
            pip.configure_nozzle_layout(
                style=PARTIAL_COLUMN,
                start= "H1",
                end = "HGFEDCBA"[channels-1]+"1"
            )
        elif channels == 8:
            pip.configure_nozzle_layout(
                style=ALL
            )
        elif channels == 1:
            pip.configure_nozzle_layout(
                style=SINGLE,
                start="H1"
            )
        if available_tips !=[]:
            if len(available_tips[0]) >= channels:
                pip.pick_up_tip(available_tips[0][channels-1])
                available_tips[0] = available_tips[0][channels:]
            else:
                remaining_tips.append(available_tips[0])
                available_tips.pop(0)
                pip.pick_up_tip(available_tips[0][channels-1])
                available_tips[0] = available_tips[0][channels:]
        
        else:
            tip = search_tip(channels)
            pip.pick_up_tip(tip)

    available_tips = partial_50.columns()
    remaining_tips = []

    ch = [5,3,2,7,6,1,8,8,6,7,4,5,3,8,8]

    for en, i in enumerate(ch):
        print(en, i)
        pick(i)
        p50_multi.drop_tip()