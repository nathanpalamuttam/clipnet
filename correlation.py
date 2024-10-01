from statistics import mean
import pandas as pd
import gzip
from scipy.stats import pearsonr

# Step 1: Loop through each seq file
for i in range(2, 8):
    csv_file = f'/home2/npp8/data/seq{i}_pause_index_run2.csv'
    bed_file = f'/fs/cbsubscb17/storage/projects/JIA_PROcap/JIA_PROcap_mapping/seq_merged/pausing_index/Seq_{i}_pausing_index.bed.gz'
    
    # Read the CSV file and extract GeneID and third column (col2 = GeneID, col3 = third column)
    csv_data = pd.read_csv(csv_file)
    print(csv_data.head())

    # Filter out rows where the third column is "undefined"
    csv_third_column = pd.to_numeric(csv_data.iloc[:, 2], errors='coerce')  # 3rd column
    csv_gene_id = csv_data.iloc[:, 1]  # 2nd column (GeneID)

    # Drop rows where the 3rd column is NaN
    csv_cleaned = pd.DataFrame({'GeneID': csv_gene_id, 'csv_col': csv_third_column}).dropna()

    print(csv_cleaned.head())

    # Step 2: Read the .bed.gz file and extract the GeneID and 7th column (col4 = GeneID, col7 = 7th column)
    bed_gene_id = []
    bed_seventh_column = []
    
    with gzip.open(bed_file, 'rt') as f:
        for line_num, line in enumerate(f, start=1):
            fields = line.strip().split()
            
            if len(fields) >= 8:
                try:
                    # Extract GeneID from col4 (3rd column in BED file) and 7th column value
                    gene_id = fields[3]  # Adjust the index based on the actual position of GeneID in the .bed file
                    bed_gene_id.append(gene_id)
                    value = float(fields[6])  # Convert 7th column to float
                    bed_seventh_column.append(value)
                except ValueError:
                    continue  # Skip lines where the conversion fails

    # Create a DataFrame for the bed data
    bed_data = pd.DataFrame({'GeneID': bed_gene_id, 'bed_col': bed_seventh_column})
    
    print(bed_data.head())

    # Step 3: Merge the two DataFrames on the 'GeneID' column
    merged_data = pd.merge(csv_cleaned, bed_data, on='GeneID')
    
    print(merged_data.head())

    # Step 4: Perform the correlation between the third column of the CSV and the 7th column of the .bed.gz file
    correlation, p_value = pearsonr(merged_data['csv_col'], merged_data['bed_col'])

    print(f'Pearson correlation: {correlation:.4f}, P-value: {p_value:.4e}')
