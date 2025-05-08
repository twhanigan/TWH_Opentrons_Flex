num_samples=10
sample_locations = []
for i in range(num_samples):
    if i < 6:  # B1 to B6
        sample_locations.append(f'B{i + 1}')
    elif i < 12:  # C1 to C6
        sample_locations.append(f'C{i - 5}')
    elif i < 18:  # D1 to D6
        sample_locations.append(f'D{i - 11}')
    else:
        break  # Stop if we exceed the number of available rows/columns
# Predefined list of letters A-H
row = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']

# Create a list of rows that repeats based on num_samples
rows = [row[i % len(row)] for i in range(num_samples)]

# Create a dynamic sample map based on the assigned sample locations
sample_map = list(map(lambda i,j :(i,j), rows, sample_locations))

filled_columns = set()
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
    filled_columns.update([1,2,3,base_column, base_column + 1, base_column + 2])
    filled_columns.update()
# After loop: compute total volume needed for all used columns
total_columns = len(filled_columns)
comb_vol = 100  # µL per well
overage_factor = 1.5
total_reagent_vol = total_columns * 8 * comb_vol * overage_factor
A_vol = total_reagent_vol * 0.50/8
B_vol = total_reagent_vol * 0.48/8
C_vol = total_reagent_vol * 0.02/8

print(A_vol,B_vol,C_vol)

# Get full columns from sample map (above)
dest_cols = [i for i in sorted(filled_columns)]
print(dest_cols)
# Dispense working reagent only to the top of filled wells


#Step 16: move plate 2 to the heater shaker and incubate at 37c
heater_shaker.open_labware_latch()
protocol.move_labware(labware=plate2, new_location=hs_adapter,use_gripper=True)
heater_shaker.close_labware_latch()
heater_shaker.set_and_wait_for_shake_speed(500)
protocol.delay(minutes=15)

#Step 17 deactivate heater shaker and temp modules
heater_shaker.deactivate_shaker()
heater_shaker.deactivate_heater()
heater_shaker.open_labware_latch()