import math
from statistics import mean
import pandas as pd
import gzip
from scipy.stats import pearsonr

# Loop through the different files (assumed 6 total in your case)
for i in range(2, 8):
    # Step 1: Read the CSV file and extract the third column, ignoring "undefined" values
    csv_file = f'/home2/npp8/data/seq{i}_pause_index_run2.csv'
    csv_data = pd.read_csv(csv_file)
    
    # Convert third column to numeric and filter out "undefined" or invalid values
    csv_third_column = pd.to_numeric(csv_data.iloc[:, 2], errors='coerce')
    csv_third_column = csv_third_column.dropna()
    
    # Step 2: Read the .bed.gz file and extract the seventh column (as this corresponds to your target)
    bed_file = f'/fs/cbsubscb17/storage/projects/JIA_PROcap/JIA_PROcap_mapping/seq_merged/pausing_index/Seq_{i}_pausing_index.bed.gz'
    bed_seventh_column = []
    
    # Open and read the .bed.gz file, extracting the seventh column (fields[6])
    with gzip.open(bed_file, 'rt') as f:
        for line_num, line in enumerate(f, start=1):
            fields = line.strip().split()
            if len(fields) >= 8:
                try:
                    value = float(fields[7])  # Extract the 7th column (0-based index = 6)
                    bed_seventh_column.append(value)
                except ValueError:
                    continue  # Skip if conversion fails (e.g., invalid data in this column)
    
    # Convert the list to a pandas Series for easier manipulation
    bed_seventh_column = pd.Series(bed_seventh_column)
    
    # print("Mean CSV Third Column:", mean(csv_third_column))
    # print("Mean BED Seventh Column:", mean(bed_seventh_column))
    
    # Combine both columns (csv_third_column and bed_seventh_column) into a DataFrame
    # Drop rows with NaN values to ensure both columns align properly
    clean_data = pd.concat([csv_third_column.reset_index(drop=True), bed_seventh_column.reset_index(drop=True)], axis=1).dropna()
    clean_data.columns = ['csv_col', 'bed_col']
    count = 0
    epsilon = 1e-5
    for index, row in clean_data.iterrows():
        
        if not math.isclose(row['csv_col'], row['bed_col'], rel_tol=epsilon, abs_tol=epsilon):
            print(f"Row {index}: csv_col = {row['csv_col']}, bed_col = {row['bed_col']}")
            if count == 5:
                break
            count += 1
            #print(f"Row {index}: csv_col = {row['csv_col']}, bed_col = {row['bed_col']}")
    #print(clean_data.head())

    # Step 3: Calculate Pearson correlation between the two columns
    if len(clean_data) > 1:  # Ensure there's enough data to calculate correlation
        correlation, p_value = pearsonr(clean_data['csv_col'], clean_data['bed_col'])
        print(f'Pearson correlation for seq{i}: {correlation:.4f}, P-value: {p_value:.4e}')
    else:
        print(f'Not enough data to calculate correlation for seq{i}')
