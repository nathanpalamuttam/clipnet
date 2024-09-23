import pickle
import csv

# Load the hashMapTSS from the pickle file
input_file = "/home2/npp8/data/hashMapTSS.pkl"
output_csv = "/home2/npp8/data/hashMapTSS_ratios.csv"

# Load the pickle object
with open(input_file, 'rb') as file:
    hashMapTSS = pickle.load(file)

# Open the output CSV file
with open(output_csv, 'w', newline='') as csvfile:
    csv_writer = csv.writer(csvfile)
    
    # Write the header to the CSV file
    csv_writer.writerow(['Key', 'GeneID', 'Pol II Ratio (TSS / TES)', 'Strand', 'TSS', 'TES'])
    
    # Iterate through the hashmap to extract required values
    for key, entries in hashMapTSS.items():
        for entry in entries:
            TSS = entry[0]
            TES = entry[1]
            pol_ii_tss = entry[2]  # # Pol II in TSS
            pol_ii_tes = entry[3]  # # Pol II in TES
            strand = entry[4]
            gene_id = entry[5]
            
            # Calculate the ratio, handling the case where TES count is 0 to avoid division by zero
            if pol_ii_tes != 0:
                denominator = (pol_ii_tes - 300) - (pol_ii_tss + 150)
                if denominator != 0:
                    ratio = (pol_ii_tss / 300) / (pol_ii_tes / denominator)
                    if ratio < 0:
                        print(pol_ii_tss)
                        print(pol_ii_tes)
                        print()
                else:
                    ratio = 'undefined'  # or some other handling for this specific case
            else:
                ratio = 'undefined' 
            # if ratio == 'undefined' or ratio == 0:
            #     print(strand)
            
            # Write the row to the CSV file
            csv_writer.writerow([key, gene_id, ratio, strand, TSS, TES])

print(f"Data successfully written to {output_csv}")
