
# Determine number of full columns to fill and convert columns to top row wells
num_full_columns = (10 + 7) // 8  # Round up to ensure all samples are covered
destination_columns = plate3.columns()[:num_full_columns]
destination_wells_col = [col[0] for col in destination_columns]  # Only use top row for multi-channel pipette
for well in destination_wells_col:
    print(well)