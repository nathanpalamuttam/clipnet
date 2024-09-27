from statistics import mean
import pandas as pd
import gzip
from scipy.stats import pearsonr

# Step 1: Read the CSV file and extract the third column, ignoring "undefined" values
csv_file = '/home2/npp8/data/seq3_pause_index_run2.csv'
csv_data = pd.read_csv(csv_file)

# Filter out rows where the third column is "undefined"
csv_third_column = pd.to_numeric(csv_data.iloc[:, 2], errors='coerce')

# Drop rows with NaN values to skip rows that had 'undefined'
csv_third_column = csv_third_column.dropna()

# Step 2: Read the .bed.gz file and extract the third column
bed_file = '/fs/cbsubscb17/storage/projects/JIA_PROcap/JIA_PROcap_mapping/seq_merged/pausing_index/Seq_3_pausing_index.bed.gz'
bed_third_column = []
count = 0
# Open and read the .bed.gz file
with gzip.open(bed_file, 'rt') as f:
    for line_num, line in enumerate(f, start=1):
        fields = line.strip().split()
        
        if len(fields) >= 8:
            try:
                count += 1
                if count <5 :
                    print(fields)
                value = float(fields[7])  # Convert the third column to float
                bed_third_column.append(value)
            except ValueError:
                continue  # Skip lines where the conversion fails
        else:
            print(f"IndexError on line {line_num}: {fields} (fields length: {len(fields)})")

print(mean(csv_third_column))
print()
print(mean(bed_third_column))
# Convert bed third column to a pandas Series
bed_third_column = pd.Series(bed_third_column)

# Step 3: Clean the data by aligning and dropping NaN values
# Combine both columns into a single DataFrame
clean_data = pd.concat([csv_third_column, bed_third_column], axis=1).dropna()
clean_data.columns = ['csv_col', 'bed_col']

# Step 4: Calculate the Pearson correlation
correlation, p_value = pearsonr(clean_data['csv_col'], clean_data['bed_col'])

print(f'Pearson correlation: {correlation:.4f}, P-value: {p_value:.4e}')
