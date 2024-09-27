import pandas as pd
import gzip
from scipy.stats import pearsonr

# Step 1: Read the CSV file and extract the third column, ignoring "undefined" values
csv_file = '/home2/npp8/data/seq3_pause_index_run2.csv'
csv_data = pd.read_csv(csv_file)

# Filter out rows where the third column is "undefined"
csv_third_column = pd.to_numeric(csv_data.iloc[:, 2], errors='coerce')  # Convert to numeric, NaNs for invalid values

# Step 2: Read the .bed.gz file and extract the third column
bed_file = '/fs/cbsubscb17/storage/projects/JIA_PROcap/JIA_PROcap_mapping/seq_merged/pausing_index/Seq_3_pausing_index.bed.gz'
bed_third_column = []
count = 0
# Open and read the .bed.gz file
with gzip.open(bed_file, 'rt') as f:
    for line in f:
        fields = line.strip().split()
        count += 1
        if count <5 :
            print(fields)
        if len(fields) >= 3:
            try:
                value = float(fields[6])  # Convert the third column to float
                bed_third_column.append(value)
            except ValueError:
                continue  # Skip lines where the conversion fails

# Convert bed third column to a pandas Series
bed_third_column = pd.Series(bed_third_column)

# Step 3: Clean the data by aligning and dropping NaN values
# Combine both columns into a single DataFrame
clean_data = pd.concat([csv_third_column, bed_third_column], axis=1).dropna()
clean_data.columns = ['csv_col', 'bed_col']

# Step 4: Calculate the Pearson correlation
correlation, p_value = pearsonr(clean_data['csv_col'], clean_data['bed_col'])

print(f'Pearson correlation: {correlation:.4f}, P-value: {p_value:.4e}')
