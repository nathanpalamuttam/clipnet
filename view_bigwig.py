from collections import defaultdict
import pickle
import pyBigWig
import gzip
# Open the BigWig file
bw = pyBigWig.open("/fs/cbsubscb17/storage/projects/JIA_PROcap/JIA_PROcap_mapping/seq_merged/Seq_dedup_QC_end_plus_merged.bw")
line_count = 0
hashMapTSS = defaultdict(list)
#this is the TSS start and end sites
with gzip.open("/fs/cbsubscb17/storage/projects/JIA_PROcap/JIA_PROcap_mapping/seq_merged/denr_greater_than_1rpb_tx.bed.gz", 'rt') as f:
    for line in f:  # adjust the range as needed to see more lines
        temp = line.strip().split("\t")
        key = temp[0]
        TSS = int(temp[1])
        TES = int(temp[2])
        strand = temp[4]
        if strand == '+':
            hashMapTSS[key].append([TSS, TES, 0, 0])
        else:
            hashMapTSS[key].append([TES, TSS, 0, 0])
print(hashMapTSS.keys())

#go through BigWig file and create hashmap where key is chromosome and value is a list
#list has 4 values: [TSS, TES, # Pol II in TSS, # Pol II in TES]

print(f"Total number of lines: {line_count}")
print("Header Information:")
print(f"Chromosomes: {bw.chroms()}")
print(f"Header Details: {bw.header()}")

# Display the first few intervals from the first chromosome
print("\nFirst few intervals:")
for chromosome in bw.chroms():
    if chromosome in hashMapTSS:
        for elem in hashMapTSS[chromosome]:
            try:
                TSS = elem[0]
                TES = elem[1]
                intervals = bw.intervals(chromosome, TSS - 150, TES - 300)
                if intervals:
                    for interval in intervals:
                        if interval[0] <= elem[0] + 150:
                            elem[2] += 1
                        else:
                            elem[3] += 1
            except Exception as e:
                print(f"Error fetching intervals: {e}")
output_file = "/home2/npp8/data/hashMapTSS.pkl"

# Save hashMapTSS to the specified file using pickle
with open(output_file, 'wb') as file:
    pickle.dump(hashMapTSS, file)

print(f"hashMapTSS has been successfully saved to {output_file}")
# Replace 'chr1' with the appropriate chromosome name present in your BigWig file
# chromosome = list(bw.chroms().keys())[0]  # Get the first chromosome in the file
# try:
#     intervals = bw.intervals(chromosome, 0, 1000)  # Fetch intervals within the first 1000 bases

#     # Print the first 10 intervals
#     for interval in intervals[:10]:  # Adjust the slice as needed
#         print(interval)
# except Exception as e:
#     print(f"Error fetching intervals: {e}")

# Close the BigWig file
bw.close()
