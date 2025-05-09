num_samples = 10
# Predefined list of letters A-H
row = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
# assign sample locations dynamically
sample_locations = []

for i in range(num_samples):
    if i < 6:  # B1 to B6
        sample_locations.append(f'B{i + 1}')
    elif i < 12:  # C1 to C6
        sample_locations.append(f'C{i - 5}')
    elif i < 18:  # D1 to D6
        sample_locations.append(f'D{i - 11}')
    elif i < 23:
        sample_locations.append(f'A{i - 16}')
    else:
        print('Too many samples')

# Generate destination wells: C1-3 → D1-3 → ... → H1-3 → A4-6 → B4-6 → ...
for index, tube in enumerate(sample_locations):
    # Row cycles C, D, E, F, G, H, A, B (skipping A/B first)
    row_index = (index + 2) % len(row)  # Start at C (index 2)
    row_letter = row[row_index]
    
    # Column block (1-3, 4-6, etc.) increments every 6 rows (C-H + A-B)
    column_block = 1 + 3 * (index // 6)  # Changes after C-H-A-B cycle
    
    destination_wells = [f'{row_letter}{column_block + i}' for i in range(3)]
    print(tube, destination_wells)
#
num_replicates = 2
loadingbuff_steps = (num_samples + 2)   # Total wells (samples + controls) * replicates
rounded = ((loadingbuff_steps) // 8 * num_replicates) + num_replicates# Round up division by 8
print(rounded)
destination_cols = [f'A{i}' for i in range(1,rounded+1)]
[print(i) for i in range(1,rounded)]
[print(destination_cols)]