# assign sample locations dynamically
num_samples = 10
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

# Iterate over the sample_map list
for index, (row, tube) in enumerate(sample_map):
    print(index,row,tube)
    if index < 8:
        base_column = 4 + (index // 8)  # This will determine the starting column for each row
    elif index< 16:
        base_column = 6 + (index // 8)
    else:
        base_column = 8 + (index //8)

    # Prepare destination wells
    destination_wells = [f'{row}{base_column + (i % 3)}' for i in range(3)]  # Generate wells like A4, A5, A6 or B4, B5, B6, etc.
    
    #Transfer the samples onto plate 2
    print(destination_wells),

