from opentrons import protocol_api
from opentrons.protocol_api import SINGLE, ALL
import pandas as pd
import numpy as np
import subprocess
from pathlib import Path
import datetime
import time
from typing import Optional, Sequence

metadata = {
    'protocolName': 'Practice Next Tip',
    'author': 'Assistant',
    'description': 'Serial dilution of BSA standard and sample processing. This includes cooling samples to 4c, heating plate to 37c with shaking and recording a video of the whole process. Place BSA Standard in A1, Lysis buffer in A2, change the number of samples and place samples in row B starting at B1. MINIMUM Sample volumen in eppendorf tubes is 40 uL. '
}

requirements = {
    "robotType": "Flex",
    "apiLevel": "2.21"
}

def next_tip(self,num_tips: int = 1,starting_tip: Optional[Well] = None,*,
nozzle_map: Optional[NozzleMapInterface] = None,) -> Optional[Well]:
    """Find the next valid well for pick-up.
    Determines the next valid start tip from which to retrieve the
    specified number of tips. There must be at least `num_tips` sequential
    wells for which all wells have tips, in the same column.
    :param num_tips: target number of sequential tips in the same column
    :type num_tips: int
    :param starting_tip: The :py:class:`.Well` from which to start search.
            for an available tip.
    :type starting_tip: :py:class:`.Well`
    :return: the :py:class:`.Well` meeting the target criteria, or None"""

    assert num_tips > 0, f"num_tips must be positive integer, but got {num_tips}"

    well_name = self._core.get_next_tip(
        num_tips=num_tips,
        starting_tip=starting_tip._core if starting_tip else None,
        nozzle_map=nozzle_map,
    )

    return self._wells_by_name[well_name] if well_name is not None else None

def run(protocol: protocol_api.ProtocolContext):
    chute = protocol.load_waste_chute()
    tips_50 = protocol.load_labware('opentrons_flex_96_filtertiprack_50ul', 'A4')
    partial_50 = protocol.load_labware(load_name="opentrons_flex_96_filtertiprack_50ul",location="A2")

    # Load pipettes
    p50_multi = protocol.load_instrument('flex_8channel_50', 'left') 
    p1000_multi = protocol.load_instrument('flex_8channel_1000', 'right')

    p50_multi.configure_nozzle_layout(style=SINGLE, start="A1",tip_racks=[partial_50])

    p50_multi.pick_up_tip()
    p50_multi.drop_tip()
    
    map2 = p50_multi.configure_nozzle_layout(style=ALL, start="A1",tip_racks=[partial_50])
    available_tip = p50_multi.next_tip(num_tips=8, starting_tip= partial_50["A1"], nozzle_map = map2)
    if available_tip is None:
        protocol.move_labware(labware=tips_50, new_location="A3", use_gripper=True)
    p50_multi.pick_up_tip()
    p50_multi.drop_tip()