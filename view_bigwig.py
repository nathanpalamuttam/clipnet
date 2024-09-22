import pyBigWig
import gzip
# Open the BigWig file
bw = pyBigWig.open("/fs/cbsubscb17/storage/projects/JIA_PROcap/JIA_PROcap_mapping/seq_merged/Seq_dedup_QC_end_plus_merged.bw")
line_count = 0
hashMapTSS = {}
#this is the TSS start and end sites
with gzip.open("/fs/cbsubscb17/storage/projects/JIA_PROcap/JIA_PROcap_mapping/seq_merged/denr_greater_than_1rpb_tx.bed.gz", 'rt') as f:
    for _ in range(5):  # adjust the range as needed to see more lines
        temp = f.readline().split(" ")
        print(temp[0])
        print(temp[1])
    for line in f:
        line_count += 1

#go through BigWig file and create hashmap where key is chromosome and value is a list
#list has 4 values: [TSS, TES, # Pol II in TSS, # Pol II in TES]

print(f"Total number of lines: {line_count}")
print("Header Information:")
print(f"Chromosomes: {bw.chroms()}")
print(f"Header Details: {bw.header()}")

# Display the first few intervals from the first chromosome
print("\nFirst few intervals:")

# Replace 'chr1' with the appropriate chromosome name present in your BigWig file
chromosome = list(bw.chroms().keys())[0]  # Get the first chromosome in the file
try:
    intervals = bw.intervals(chromosome, 0, 1000)  # Fetch intervals within the first 1000 bases

    # Print the first 10 intervals
    for interval in intervals[:10]:  # Adjust the slice as needed
        print(interval)
except Exception as e:
    print(f"Error fetching intervals: {e}")

# Close the BigWig file
bw.close()
