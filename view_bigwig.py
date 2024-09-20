import pyBigWig

# Open the BigWig file
bw = pyBigWig.open("Seq_dedup_QC_end_minus_merged.bw")

# Print the header information
print("Header Information:")
print(f"Chromosomes: {bw.chroms()}")
print(f"Total Items: {bw.header()['nBasesCovered']}")

# Print the first few entries from the first chromosome
print("\nFirst few entries:")

# Replace 'chr1' with the appropriate chromosome name present in your BigWig file
chromosome = list(bw.chroms().keys())[0]  # Get the first chromosome
entries = bw.entries(chromosome, 0, 1000)  # Read first 1000 bases of the chromosome

# Display the first few entries
for entry in entries[:10]:  # Adjust number to display more or fewer entries
    print(entry)

# Close the file
bw.close()