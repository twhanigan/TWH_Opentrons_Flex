from opentrons import protocol_api
from opentrons.protocol_api import SINGLE, PARTIAL_COLUMN, ALL
from opentrons import types


metadata = {
    'protocolName': 'Partial Tip Racks',
    'author': 'Anurag Kanase <anurag.kanase@opentrons.com>',
    'description': 'A protocol for using partial tip racks.'
}

requirements = {'apiLevel': '2.21', 'robotType': 'Flex'}

def add_parameters(parameters: protocol_api.Parameters): 

    parameters.add_str(
        display_name="First Tip row",
        variable_name="first_tip_row",
        default="A",
        choices=[{"display_name": row, "value": row} for row in "ABCDEFGH"],
        description="Select first tip row to pickup from."
    )
    parameters.add_int(
        display_name="First Tip column",
        variable_name="first_tip_column",
        default=1,minimum=1,maximum=96,
        description="Select first tip column to pickup from.")
    


def run(protocol: protocol_api.ProtocolContext):

    trash = protocol.load_trash_bin("A3")

    first_tip_column = protocol.params.first_tip_column
    first_tip_row = protocol.params.first_tip_row
    start_tip = first_tip_row + str(first_tip_column)

    tiprack = protocol.load_labware("opentrons_flex_96_tiprack_1000ul", "B1")

    m1000 = protocol.load_instrument("flex_8channel_1000", "left")

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



    def pick(channels, pip = m1000):
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

    
    available_tips = tiprack.columns()
    remaining_tips = []

    ch = [5,3,2,7,6,1,8,8,6,7,4,5,3,8,8]

    for en, i in enumerate(ch):
        print(en, i)
        pick(i)
        m1000.drop_tip()

    # try:
    #     for en, i in enumerate(ch):
    #         print(en, i)
    #         pick(i)
    #         m1000.drop_tip()
    # except:
    #     print("Available tips:")
    #     print(available_tips)
    #     print("Remaining tips:")
    #     print(remaining_tips)



    
    
    