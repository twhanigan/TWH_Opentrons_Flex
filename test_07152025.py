# assign sample locations dynamically

# Enter the number of samples 
speed= 0.5
target_concentration = 1
num_samples = 15
num_replicates = 3
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
numtotalSamples = num_samples + 2
reaction_vol = 25
mastermix_vol = 12.5
primer_vol = 1
sample_vol = 2.5
water_vol = 9

#[print(well) for well in sample_well_map["positive_control"]]
# Transfer positive control to A1
#print("negative control")
#[print(well) for well in sample_well_map["neg_control"]]
#print("samples")
# Transfer samples
for sample_idx in range(num_samples):
    tube_loc = sample_locations[sample_idx]
    wells = sample_well_map[f"sample_{sample_idx+1}"]
    #[print(well) for well in wells]
    print(tube_loc)