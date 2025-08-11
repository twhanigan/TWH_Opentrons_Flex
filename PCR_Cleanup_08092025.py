from opentrons import protocol_api
from opentrons.protocol_api import SINGLE, ALL
import subprocess

metadata = {
    'protocolName': 'PCR cleanup for NGS',
    'author': 'Assistant',
    'description': 'PCR Amplification of OXA1L'
}

requirements = {
    "robotType": "Flex",
    "apiLevel": "2.21"
}

import math

def add_parameters(parameters):
    # ======================== RUNTIME PARAMETERS ========================
    parameters.add_int(
        display_name="Number of samples",
        variable_name="NUM_SAMPLES",
        default=24, minimum=1, maximum=96,
        description="Total number of samples to process (fills by column with an 8‑channel)."
    )
    parameters.add_int(
        display_name="Input Volume (µl)",
        variable_name="INPUTVOLUME",
        default=20, minimum=10, maximum=100,
        description="Input volume."
    )
    parameters.add_float(
        display_name="Bead ratio (x)",
        variable_name="BEADRATIO",
        default=2, minimum=0.5, maximum=2,
        description='Bead Ratio relative to input volume in "x".'
    )

# ---------- helpers ----------
ROWS = ['A','B','C','D','E','F','G','H']

def wells_for_n(n, start_col=1):
    """Return first n wells, filling by columns from start_col (A..H then next col)."""
    out = []
    col = start_col
    while len(out) < n:
        for r in ROWS:
            if len(out) == n:
                break
            out.append(f"{r}{col}")
        col += 1
    return out

def split_by_column(well_names):
    """Map column index -> [wells in that column], preserving A..H order."""
    bycol = {}
    for w in well_names:
        # well like 'A12' -> 12
        col = int(w[1:])
        bycol.setdefault(col, []).append(w)
    # ensure A..H order per column
    for c in bycol:
        bycol[c] = [f"{r}{c}" for r in ROWS if f"{r}{c}" in bycol[c]]
    return bycol

def run(protocol: protocol_api.ProtocolContext):
    # ======== PULL PARAMS ========
    # (Assuming these are provided by Opentrons param system)
    # NUM_SAMPLES, INPUTVOLUME, BEADRATIO, BEADRATIO_2, RES_TYPE_96x, etc. exist.

    # ======== DERIVED FROM NUM_SAMPLES ========
    COLUMNS = math.ceil(NUM_SAMPLES / 8)   # <- replaces user-entered COLUMNS everywhere
    samples_src = wells_for_n(NUM_SAMPLES, start_col=1)  # source wells (plate 1, cols 1..)
    # A "second set" of wells (same count) immediately after sources, used for double-sided &/or final placement
    second_set_start = COLUMNS + 1
    samples_next = wells_for_n(NUM_SAMPLES, start_col=second_set_start)

    # If you need explicit per-step lists like in your code:
    column_1_list = samples_src                    # original sample positions
    column_2_list = samples_next                   # second pass / second column set
    # final elution targets (if NEW_PLATE: plate 2 col1.., else plate 1 after sources)
    if NEW_PLATE:
        column_3_list = wells_for_n(NUM_SAMPLES, start_col=1)  # on sample_plate_2
    else:
        column_3_list = wells_for_n(NUM_SAMPLES, start_col=second_set_start)  # on sample_plate_1

    # (Optional) access by actual column if you need to iterate column-wise:
    src_by_col = split_by_column(column_1_list)
    next_by_col = split_by_column(column_2_list)

    # ======== LIQUID ESTIMATES ========
    Sample_Volume = INPUTVOLUME
    AMPure_Volume = COLUMNS * ((BEADRATIO + BEADRATIO_2) * INPUTVOLUME)
    ETOH_Volume   = 8 * COLUMNS * (150 * 2)
    RSB_Volume    = 8 * COLUMNS * (RESUSPENSION)

    TotalRowLetters = ROWS[:]  # ['A'..'H']

    # ======== DEFINING LIQUIDS (unchanged) ========
    AMPure = protocol.define_liquid(name="Beads", description="AMPure Beads", display_color="#704848")
    EtOH = protocol.define_liquid(name="EtOH", description="80% Ethanol", display_color="#9ACECB")
    RSB = protocol.define_liquid(name="RSB", description="Resuspension Buffer", display_color="#00FFF2")
    Liquid_trash_well = protocol.define_liquid(name="Liquid_trash_well", description="Liquid Trash", display_color="#9B9B9B")
    Sample = protocol.define_liquid(name="Sample", description="Sample", display_color="#52AAFF")
    Final_Sample = protocol.define_liquid(name="Final_Sample", description="Final Sample", display_color="#82A9CF")

    # ======== RESERVOIR LOADING (uses derived COLUMNS exactly like before) ========
    if RES_TYPE_96x == False:
        reservoir.wells_by_name()['A1'].load_liquid(liquid=AMPure, volume=AMPure_Volume*8+1000)
        if ETOH_Volume <= 15000:
            reservoir.wells_by_name()['A4'].load_liquid(liquid=EtOH, volume=15000)
        if 15000 <= ETOH_Volume < 30000:
            reservoir.wells_by_name()['A4'].load_liquid(liquid=EtOH, volume=15000)
            reservoir.wells_by_name()['A5'].load_liquid(liquid=EtOH, volume=15000)
        if ETOH_Volume >= 30000:
            reservoir.wells_by_name()['A4'].load_liquid(liquid=EtOH, volume=15000)
            reservoir.wells_by_name()['A5'].load_liquid(liquid=EtOH, volume=15000)
            reservoir.wells_by_name()['A6'].load_liquid(liquid=EtOH, volume=15000)
        reservoir.wells_by_name()['A7'].load_liquid(liquid=RSB, volume=RSB_Volume*8)
        reservoir.wells_by_name()['A10'].load_liquid(liquid=Liquid_trash_well, volume=0)
        reservoir.wells_by_name()['A11'].load_liquid(liquid=Liquid_trash_well, volume=0)
        reservoir.wells_by_name()['A12'].load_liquid(liquid=Liquid_trash_well, volume=0)
    else:
        # if you previously looped over 8 row-letters, you can still do that,
        # but only as many columns as COLUMNS (derived from NUM_SAMPLES)
        UsedColumnLetters = ROWS[:]  # unchanged
        for loop, X in enumerate(UsedColumnLetters):
            reservoir.wells_by_name()[X+'1'].load_liquid(liquid=AMPure, volume=AMPure_Volume)
            if ETOH_Volume >= 15000:
                reservoir.wells_by_name()[X+'4'].load_liquid(liquid=EtOh, volume=15000)
            if 15000 <= ETOH_Volume < 30000:
                reservoir.wells_by_name()[X+'4'].load_liquid(liquid=EtOH, volume=15000)
                reservoir.wells_by_name()[X+'5'].load_liquid(liquid=EtOH, volume=15000)
            if ETOH_Volume >= 30000:
                reservoir.wells_by_name()[X+'4'].load_liquid(liquid=EtOH, volume=15000)
                reservoir.wells_by_name()[X+'5'].load_liquid(liquid=EtOH, volume=15000)
                reservoir.wells_by_name()[X+'6'].load_liquid(liquid=EtOH, volume=15000)
            reservoir.wells_by_name()[X+'7'].load_liquid(liquid=RSB, volume=RSB_Volume)
            reservoir.wells_by_name()[X+'10'].load_liquid(liquid=Liquid_trash_well, volume=0)
            reservoir.wells_by_name()[X+'11'].load_liquid(liquid=Liquid_trash_well, volume=0)
            reservoir.wells_by_name()[X+'12'].load_liquid(liquid=Liquid_trash_well, volume=0)

    # ======== SAMPLE/DEST PLATE LIQUID "LOAD" MARKERS (now respect NUM_SAMPLES) ========
    # Mark sources
    for w in column_1_list:
        sample_plate_1.wells_by_name()[w].load_liquid(liquid=Sample, volume=Sample_Volume)

    # Mark destinations on same plate AFTER sources (if not NEW_PLATE)
    if not NEW_PLATE:
        for w in column_3_list:
            sample_plate_1.wells_by_name()[w].load_liquid(liquid=Final_Sample, volume=0)
    else:
        for w in column_3_list:
            sample_plate_2.wells_by_name()[w].load_liquid(liquid=Final_Sample, volume=0)

    # =================================================================================================
    # ========================================= PROTOCOL START ========================================
    # =================================================================================================
    if ONDECK_THERMO == True: thermocycler.open_lid()
    if ONDECK_HEATERSHAKER == True: heatershaker.open_labware_latch()
    protocol.pause("Ready")
    if ONDECK_HEATERSHAKER == True: heatershaker.close_labware_latch() 
    Liquid_trash = Liquid_trash_well_1
    # =================================================================================================
    # ========================================= PROTOCOL START ========================================
    # =================================================================================================
    
    if STEP_CLEANUP == 1:
        protocol.comment('==============================================')
        protocol.comment('--> Cleanup 1')
        protocol.comment('==============================================')

        protocol.comment('--> Transferring Sample to CleanupPlate')
        TransferSup = int(INPUTVOLUME/math.ceil(INPUTVOLUME/50))
        p50.flow_rate.aspirate = p50_flow_rate_aspirate_default*0.25
        p50.flow_rate.dispense = p50_flow_rate_dispense_default*0.75
        p50.flow_rate.blow_out = p50_flow_rate_blow_out_default*0.5
        #=============================
        for loop, X in enumerate(column_1_list):
            TipCheck(50,None,loop,None)
            for Y in range(math.ceil(INPUTVOLUME/50)):
                p50.aspirate(TransferSup, sample_plate_1[X].bottom(z=PCRPlate_Z_offset+1))
                p50.dispense(TransferSup, CleanupPlate[column_1_list[loop]].bottom(z=Deepwell_Z_offset+1))
            TipDone(50,None,loop,None)
        #=============================

        protocol.comment('--> ADDING AMPure')
        SampleVol = INPUTVOLUME
        AMPureMixRPM = SHAKE_RPM
        AMPureMixTime = SHAKE_TIME*60
        AMPurePremix = 3
        p1000.flow_rate.aspirate = p1000_flow_rate_aspirate_default*0.5
        p1000.flow_rate.dispense = p1000_flow_rate_dispense_default*0.5
        p1000.flow_rate.blow_out = p1000_flow_rate_blow_out_default*0.5
        #=============================
        for loop, X in enumerate(column_1_list):
            AMPureVol = (BEADRATIO*INPUTVOLUME)
            TipCheck(200,None,loop,None)
            p1000.aspirate(AMPureVol, AMPure.bottom(z=Deepwell_Z_offset+1.5), rate=0.25)
            p1000.dispense(AMPureVol, AMPure.bottom(z=Deepwell_Z_offset+1.5), rate=0.25)
            p1000.aspirate(AMPureVol+3 if INPUTVOLUME < 100 else AMPureVol, AMPure.bottom(z=Deepwell_Z_offset+1), rate=0.25)
            if INPUTVOLUME < 100:
                p1000.dispense(3, AMPure.bottom(z=Deepwell_Z_offset+1), rate=0.25)
            p1000.default_speed = 5
            p1000.move_to(AMPure.top(z=-3))
            #=====Reservoir Tip Touch========
            p1000.default_speed = 100
            p1000.move_to(AMPure.top().move(types.Point(x=4,z=-3)))
            p1000.move_to(AMPure.top().move(types.Point(x=-4,z=-3)))
            p1000.default_speed = 400
            #================================  
            p1000.dispense(AMPureVol, CleanupPlate[X].bottom(z=Deepwell_Z_offset+3), rate=0.5)
            p1000.move_to(CleanupPlate[X].bottom(z=Deepwell_Z_offset+7))
            p1000.default_speed = 100
            p1000.move_to(CleanupPlate[X].bottom(z=Deepwell_Z_offset+1))
            if TIP_MIX == True:
                AmpureMixRep = 10
            if TIP_MIX == False:
                AmpureMixRep = 2
            for Mix in range(AmpureMixRep):
                p1000.aspirate((SampleVol+AMPureVol)-10, rate=0.5)
                protocol.delay(seconds=1)
                p1000.dispense((SampleVol+AMPureVol)-10, rate=0.5)
                Mix += 1
            p1000.move_to(CleanupPlate[X].bottom(z=Deepwell_Z_offset+7))
            p1000.default_speed = 400
            p1000.move_to(CleanupPlate[X].top(z=-3))
            protocol.delay(seconds=1)
            p1000.blow_out(CleanupPlate[X].top(z=-3))
            p1000.touch_tip(speed=100)
            p1000.move_to(CleanupPlate[X].top(z=5))
            p1000.move_to(CleanupPlate[X].top(z=0))
            p1000.move_to(CleanupPlate[X].top(z=5))
            TipDone(200,None,loop,None)
        #=============================
        if TIP_MIX == False:
            heatershaker.set_and_wait_for_shake_speed(rpm=AMPureMixRPM)
            protocol.delay(AMPureMixTime)
            heatershaker.deactivate_shaker()

        #============================================================================================
        # GRIPPER MOVE CleanupPlate FROM HEATER SHAKER TO MAG PLATE
        if ONDECK_HEATERSHAKER == True:
            heatershaker.open_labware_latch()
            protocol.move_labware(labware=CleanupPlate,new_location=mag_block,use_gripper=USE_GRIPPER,pick_up_offset=hs_pick_up_offset,drop_offset=mb_drop_offset)
            heatershaker.close_labware_latch()
        else:
            protocol.move_labware(labware=CleanupPlate,new_location=mag_block,use_gripper=USE_GRIPPER,pick_up_offset=deck_pick_up_offset,drop_offset=mb_drop_offset)
        #============================================================================================

        protocol.delay(minutes=4)

        if DOUBLEDSIDED == True:
            protocol.comment('--> Transferring Sample to CleanupPlate')
            #TransferReps = math.ceil((INPUTVOLUME+(BEADRATIO*INPUTVOLUME))/50)
            #TransferSup = 50 if int((INPUTVOLUME+(BEADRATIO*INPUTVOLUME))/TransferReps)+3 >=50 else int((INPUTVOLUME+(BEADRATIO*INPUTVOLUME))/TransferReps)+3
            p50.flow_rate.aspirate = p50_flow_rate_aspirate_default*0.25
            p50.flow_rate.dispense = p50_flow_rate_dispense_default*0.75
            p50.flow_rate.blow_out = p50_flow_rate_blow_out_default*0.5
            #=============================
            for loop, X in enumerate(column_1_list):
                TransferReps = math.ceil((INPUTVOLUME+(BEADRATIO*INPUTVOLUME))/50)
                TransferSup = 50 if int((INPUTVOLUME+(BEADRATIO*INPUTVOLUME))/TransferReps)+3 >=50 else int((INPUTVOLUME+(BEADRATIO*INPUTVOLUME))/TransferReps)+3
                TipCheck(50,None,loop,None)
                for Y in range(TransferReps):
                    p50.aspirate(TransferSup, CleanupPlate[X].bottom(z=Deepwell_Z_offset+1))
                    p50.dispense(TransferSup, CleanupPlate[column_2_list[loop]].bottom(z=Deepwell_Z_offset+1.5))
                TipDone(50,None,loop,None)
            #=============================

            #============================================================================================
            # GRIPPER MOVE CleanupPlate FROM MAG PLATE TO HEATER SHAKER
            if ONDECK_HEATERSHAKER == True:
                heatershaker.open_labware_latch()
                protocol.move_labware(labware=CleanupPlate,new_location=hs_adapter,use_gripper=USE_GRIPPER,pick_up_offset=mb_pick_up_offset,drop_offset=hs_drop_offset)
                heatershaker.close_labware_latch()
            else:
                protocol.move_labware(labware=CleanupPlate,new_location='D1',use_gripper=USE_GRIPPER,pick_up_offset=deck_pick_up_offset,drop_offset=deck_drop_offset)
            #============================================================================================

            protocol.comment('--> ADDING AMPure (0.8x)')
            AMPureVol = (BEADRATIO_2*INPUTVOLUME)
            SampleVol = INPUTVOLUME
            AMPureMixRPM = SHAKE_RPM
            AMPureMixTime = SHAKE_TIME*60
            AMPurePremix = 3
            p1000.flow_rate.aspirate = p1000_flow_rate_aspirate_default*0.5
            p1000.flow_rate.dispense = p1000_flow_rate_dispense_default*0.5
            p1000.flow_rate.blow_out = p1000_flow_rate_blow_out_default*0.5
            #=============================
            for loop, X in enumerate(column_2_list):
                TipCheck(200,None,loop,None)
                p1000.aspirate(AMPureVol, AMPure.bottom(z=Deepwell_Z_offset+1.5), rate=0.25)
                p1000.dispense(AMPureVol, AMPure.bottom(z=Deepwell_Z_offset+1.5), rate=0.25)
                p1000.aspirate(AMPureVol+3, AMPure.bottom(z=Deepwell_Z_offset+1.5), rate=0.25)
                p1000.dispense(3, AMPure.bottom(z=Deepwell_Z_offset+1.5), rate=0.25)
                p1000.default_speed = 5
                p1000.move_to(AMPure.top(z=-3))
                #=====Reservoir Tip Touch========
                p1000.default_speed = 100
                p1000.move_to(AMPure.top().move(types.Point(x=4,z=-3)))
                p1000.move_to(AMPure.top().move(types.Point(x=-4,z=-3)))
                p1000.default_speed = 400
                #================================  
                p1000.dispense(AMPureVol, CleanupPlate[X].bottom(z=Deepwell_Z_offset+3), rate=0.5)
                p1000.move_to(CleanupPlate[X].bottom(z=Deepwell_Z_offset+7))
                p1000.default_speed = 100
                p1000.move_to(CleanupPlate[X].bottom(z=Deepwell_Z_offset+3))
                if TIP_MIX == True:
                    AmpureMixRep = 10
                if TIP_MIX == False:
                    AmpureMixRep = 2
                for Mix in range(AmpureMixRep):
                    p1000.aspirate(70, rate=0.5)
                    p1000.move_to(CleanupPlate[X].bottom(z=Deepwell_Z_offset+1))
                    p1000.aspirate(20, rate=0.5)
                    p1000.dispense(20, rate=0.5)
                    p1000.move_to(CleanupPlate[X].bottom(z=Deepwell_Z_offset+3))
                    p1000.dispense(70, rate=0.5)
                    Mix += 1
                p1000.move_to(CleanupPlate[X].bottom(z=Deepwell_Z_offset+7))
                p1000.default_speed = 400
                p1000.move_to(CleanupPlate[X].top(z=-3))
                protocol.delay(seconds=1)
                p1000.blow_out(CleanupPlate[X].top(z=-3))
                p1000.touch_tip(speed=100)
                p1000.move_to(CleanupPlate[X].top(z=5))
                p1000.move_to(CleanupPlate[X].top(z=0))
                p1000.move_to(CleanupPlate[X].top(z=5))
                TipDone(200,None,loop,None)
            #=============================
            if TIP_MIX == False:
                heatershaker.set_and_wait_for_shake_speed(rpm=AMPureMixRPM)
                protocol.delay(AMPureMixTime)
                heatershaker.deactivate_shaker()

            #============================================================================================
            # GRIPPER MOVE CleanupPlate FROM HEATER SHAKER TO MAG PLATE
            if ONDECK_HEATERSHAKER == True:
                heatershaker.open_labware_latch()
                protocol.move_labware(labware=CleanupPlate,new_location=mag_block,use_gripper=USE_GRIPPER,pick_up_offset=hs_pick_up_offset,drop_offset=mb_drop_offset)
                heatershaker.close_labware_latch()
            else:
                protocol.move_labware(labware=CleanupPlate,new_location=mag_block,use_gripper=USE_GRIPPER,pick_up_offset=deck_pick_up_offset,drop_offset=deck_drop_offset)
            #============================================================================================

        protocol.comment('--> Removing Supernatant')
        RemoveSup = 200
        ActualRemoveSup = 8*(BEADRATIO*INPUTVOLUME) + INPUTVOLUME
        p1000.flow_rate.aspirate = p1000_flow_rate_aspirate_default*0.5
        p1000.flow_rate.dispense = p1000_flow_rate_dispense_default*0.5
        p1000.flow_rate.blow_out = p1000_flow_rate_blow_out_default*0.5
        #=============================
        for loop, X in enumerate(column_2_list):
            TipCheck(200,'REUSE',loop,'ETOH')
            p1000.move_to(CleanupPlate[X].bottom(z=Deepwell_Z_offset+2))
            p1000.aspirate(RemoveSup-100)
            protocol.delay(minutes=0.1)
            p1000.move_to(CleanupPlate[X].bottom(z=Deepwell_Z_offset+1))
            p1000.aspirate(100)
            p1000.default_speed = 200
            p1000.move_to(CleanupPlate[X].top(z=2))
            #======L Waste Volume Check======
            WASTEVOL+=(ActualRemoveSup*8)
            protocol.comment('Adding '+str((ActualRemoveSup*8))+'ul tp '+str(WASTEVOL))
            if WASTEVOL <14400:
                Liquid_trash = Liquid_trash_well_1
            if WASTEVOL >=14400 and WASTEVOL <28800:
                Liquid_trash = Liquid_trash_well_2
            if WASTEVOL >=28800:
                Liquid_trash = Liquid_trash_well_3
            #================================ 
            p1000.dispense(200, Liquid_trash.top(z=0))
            protocol.delay(minutes=0.1)
            p1000.blow_out()
            p1000.default_speed = 400
            p1000.move_to(Liquid_trash.top(z=-5))
            p1000.move_to(Liquid_trash.top(z=0))
            TipDone(200,'REUSE',loop,'ETOH')
        #=============================

        for X in range(2):
            protocol.comment('--> ETOH Wash')
            ETOHMaxVol = 200
            p1000.flow_rate.aspirate = p1000_flow_rate_aspirate_default*0.5
            p1000.flow_rate.dispense = p1000_flow_rate_dispense_default
            p1000.flow_rate.blow_out = p1000_flow_rate_blow_out_default
            #=============================
            if ETOH_1_AirMultiDis == True:
                TipCheck(200,None,loop,None)
                for loop, X in enumerate(column_2_list):
                    p1000.aspirate(ETOHMaxVol, reservoir.wells_by_name()[ETOH_list[loop]].bottom(z=Deepwell_Z_offset+1.5))
                    p1000.move_to(reservoir.wells_by_name()[ETOH_list[loop]].top(z=0))
                    p1000.move_to(reservoir.wells_by_name()[ETOH_list[loop]].top(z=-5))
                    p1000.move_to(reservoir.wells_by_name()[ETOH_list[loop]].top(z=0))
                    #=====Reservoir Tip Touch========
                    p1000.default_speed = 100
                    p1000.move_to(reservoir.wells_by_name()[ETOH_list[loop]].top().move(types.Point(x=4,z=-3)))
                    p1000.move_to(reservoir.wells_by_name()[ETOH_list[loop]].top().move(types.Point(x=-4,z=-3)))
                    p1000.default_speed = 400
                    #================================ 
                    p1000.move_to(CleanupPlate[X].top(z=7))
                    p1000.dispense(ETOHMaxVol)
                    protocol.delay(minutes=0.1)
                    p1000.blow_out(CleanupPlate[X].top(z=2))
                    p1000.move_to(CleanupPlate[X].top(z=5))
                    p1000.move_to(CleanupPlate[X].top(z=2))
                    p1000.move_to(CleanupPlate[X].top(z=5))
                TipDone(200,None,loop,None)
            else:
                for loop, X in enumerate(column_2_list):
                    TipCheck(200,None,loop,None)
                    p1000.aspirate(ETOHMaxVol, reservoir.wells_by_name()[ETOH_list[loop]].bottom(z=Deepwell_Z_offset+1.5))
                    p1000.move_to(reservoir.wells_by_name()[ETOH_list[loop]].top(z=0))
                    p1000.move_to(reservoir.wells_by_name()[ETOH_list[loop]].top(z=-5))
                    p1000.move_to(reservoir.wells_by_name()[ETOH_list[loop]].top(z=0))
                    #=====Reservoir Tip Touch========
                    p1000.default_speed = 100
                    p1000.move_to(reservoir.wells_by_name()[ETOH_list[loop]].top().move(types.Point(x=4,z=-3)))
                    p1000.move_to(reservoir.wells_by_name()[ETOH_list[loop]].top().move(types.Point(x=-4,z=-3)))
                    p1000.default_speed = 400
                    #================================ 
                    p1000.move_to(CleanupPlate[X].top(z=-5))
                    p1000.dispense(ETOHMaxVol)
                    protocol.delay(minutes=0.1)
                    p1000.blow_out()
                    p1000.move_to(CleanupPlate[X].top(z=5))
                    p1000.move_to(CleanupPlate[X].top(z=0))
                    p1000.move_to(CleanupPlate[X].top(z=5))
                    TipDone(200,None,loop,None)
            #=============================
            protocol.delay(minutes=0.5)
            
            protocol.comment('--> Remove ETOH Wash')
            RemoveSup = 200
            ActualRemoveSup = 200
            p1000.flow_rate.aspirate = p1000_flow_rate_aspirate_default*0.5
            p1000.flow_rate.dispense = p1000_flow_rate_dispense_default
            p1000.flow_rate.blow_out = p1000_flow_rate_blow_out_default
            #=============================
            for loop, X in enumerate(column_2_list):
                TipCheck(200,'REUSE',loop,'ETOH')
                p1000.move_to(CleanupPlate[X].bottom(z=Deepwell_Z_offset+2))
                p1000.aspirate(RemoveSup-100)
                protocol.delay(minutes=0.1)
                p1000.move_to(CleanupPlate[X].bottom(z=Deepwell_Z_offset+1))
                p1000.aspirate(100)
                p1000.default_speed = 5
                p1000.move_to(CleanupPlate[X].bottom(z=Deepwell_Z_offset+7))
                p1000.default_speed = 200
                #======L Waste Volume Check======
                WASTEVOL+=(ActualRemoveSup*8)
                protocol.comment('Adding '+str((ActualRemoveSup*8))+'ul tp '+str(WASTEVOL))
                if WASTEVOL <14400:
                    Liquid_trash = Liquid_trash_well_1
                if WASTEVOL >=14400 and WASTEVOL <28800:
                    Liquid_trash = Liquid_trash_well_2
                if WASTEVOL >=28800:
                    Liquid_trash = Liquid_trash_well_3
                #================================ 
                p1000.dispense(200, Liquid_trash.top(z=0))
                protocol.delay(minutes=0.1)
                p1000.blow_out()
                p1000.default_speed = 400
                p1000.move_to(Liquid_trash.top(z=-5))
                p1000.move_to(Liquid_trash.top(z=0))
                TipDone(200,'REUSE',loop,'ETOH')
            #=============================

        protocol.delay(minutes=ETOH_TIME)

        for batch in range(BATCHREP):
            protocol.comment('--> Removing Residual Wash')
            p50.flow_rate.aspirate = p50_flow_rate_aspirate_default*0.5
            p50.flow_rate.dispense = p50_flow_rate_dispense_default*0.5
            p50.flow_rate.blow_out = p50_flow_rate_blow_out_default*0.5
            #=============================
            for loop, X in enumerate(column_2_list[batch*(int(len(column_2_list)/BATCHREP)):(1+batch)*int(len(column_2_list)/BATCHREP)]):
                TipCheck(50,None,loop,None)
                p50.move_to(CleanupPlate[X].bottom(z=Deepwell_Z_offset+1))
                p50.aspirate(50, rate=0.25)
                TipDone(50,None,loop,None)
            #=============================

            protocol.delay(minutes=0.5)
            
            #============================================================================================
            # GRIPPER MOVE CleanupPlate FROM MAG PLATE TO HEATER SHAKER
            if ONDECK_HEATERSHAKER == True:
                heatershaker.open_labware_latch()
                protocol.move_labware(labware=CleanupPlate,new_location=hs_adapter,use_gripper=USE_GRIPPER,pick_up_offset=mb_pick_up_offset,drop_offset=hs_drop_offset)
                heatershaker.close_labware_latch()
            else:
                protocol.move_labware(labware=CleanupPlate,new_location='D1',use_gripper=USE_GRIPPER,pick_up_offset=deck_pick_up_offset,drop_offset=deck_drop_offset)
            #============================================================================================
            
            protocol.comment('--> Adding RSB')
            RSBVol = RESUSPENSION+2
            RSBMixRPM = SHAKE_RPM
            RSBMixTime = SHAKE_TIME*60
            p50.flow_rate.aspirate = p50_flow_rate_aspirate_default*0.5
            p50.flow_rate.dispense = p50_flow_rate_dispense_default*0.5
            p50.flow_rate.blow_out = p50_flow_rate_blow_out_default*0.5
            #=============================
            for loop, X in enumerate(column_2_list[batch*(int(len(column_2_list)/BATCHREP)):(1+batch)*int(len(column_2_list)/BATCHREP)]):
                TipCheck(50,'REUSE',loop,None)
                p50.aspirate(RSBVol, RSB.bottom(z=Deepwell_Z_offset+1.5))
                p50.move_to(CleanupPlate.wells_by_name()[X].bottom(z=Deepwell_Z_offset+2))
                p50.dispense(RSBVol)
                if TIP_MIX == True:
                    RSBMixRep = 10
                if TIP_MIX == False:
                    RSBMixRep = 2
                p50.move_to(CleanupPlate[X].bottom(z=Deepwell_Z_offset+1))
                for Mix in range(RSBMixRep):
                    p50.aspirate(RSBVol-5, rate=0.5)
                    p50.dispense(RSBVol-5, rate=0.5)
                    Mix += 1
                p50.blow_out(CleanupPlate.wells_by_name()[X].top(z=3))
                TipDone(50,'REUSE',loop,None)
            #=============================
            if TIP_MIX == False:
                heatershaker.set_and_wait_for_shake_speed(rpm=RSBMixRPM)
                protocol.delay(RSBMixTime)
                heatershaker.deactivate_shaker()
            
            #============================================================================================
            # GRIPPER MOVE CleanupPlate FROM HEATER SHAKER TO MAG PLATE
            if ONDECK_HEATERSHAKER == True:
                heatershaker.open_labware_latch()
                protocol.move_labware(labware=CleanupPlate,new_location=mag_block,use_gripper=USE_GRIPPER,pick_up_offset=hs_pick_up_offset,drop_offset=mb_drop_offset)
                heatershaker.close_labware_latch()
            else:
                protocol.move_labware(labware=CleanupPlate,new_location=mag_block,use_gripper=USE_GRIPPER,pick_up_offset=deck_pick_up_offset,drop_offset=deck_drop_offset)
            #============================================================================================
            
        protocol.delay(minutes=3)

        protocol.comment('--> Transferring Supernatant')
        TransferSup = RESUSPENSION
        p50.flow_rate.aspirate = p50_flow_rate_aspirate_default*0.5
        p50.flow_rate.dispense = p50_flow_rate_dispense_default*0.5
        p50.flow_rate.blow_out = p50_flow_rate_blow_out_default*0.5
        #=============================
        for loop, X in enumerate(column_2_list):
            TipCheck(50,'REUSE',loop,None)
            p50.move_to(CleanupPlate[X].bottom(z=Deepwell_Z_offset+1))
            p50.aspirate(TransferSup)
            if NEW_PLATE == True:
                p50.dispense(TransferSup, sample_plate_2[column_3_list[loop]].bottom(z=PCRPlate_Z_offset+1))
                p50.blow_out(sample_plate_2[column_3_list[loop]].top(z=PCRPlate_Z_offset-3))
            else:
                p50.dispense(TransferSup, sample_plate_1[column_3_list[loop]].bottom(z=PCRPlate_Z_offset+1))
                p50.blow_out(sample_plate_1[column_3_list[loop]].top(z=PCRPlate_Z_offset-3))
            TipDone(50,'REUSE',loop,None)
        #=============================


    # =================================================================================================
    # ========================================== PROTOCOL END =========================================
    # =================================================================================================
    if ONDECK_HEATERSHAKER == True: heatershaker.open_labware_latch()