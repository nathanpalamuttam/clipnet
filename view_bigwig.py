import pyBigWig

# Open the BigWig file
bw = pyBigWig.open("/fs/cbsubscb17/storage/projects/JIA_PROcap/JIA_PROcap_mapping/seq_merged/Seq_dedup_QC_end_plus_merged.bw")

# Print header information
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
