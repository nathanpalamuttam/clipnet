from collections import defaultdict
import pickle
import pyBigWig
import gzip
import pickle
import csv
# Open the BigWig file
for i in range(2, 8):
    bw = pyBigWig.open(f"/fs/cbsubscb17/storage/projects/JIA_PROcap/JIA_PROcap_mapping/seq_merged/Seq_{i}/Seq_{i}_dedup_QC_end_plus.bw")
    bw_neg = pyBigWig.open(f"/fs/cbsubscb17/storage/projects/JIA_PROcap/JIA_PROcap_mapping/seq_merged/Seq_{i}/Seq_{i}_dedup_QC_end_minus.bw")

    line_count = 0
    hashMapTSS = defaultdict(list)
    #this is the TSS start and end sites
    with gzip.open("/fs/cbsubscb17/storage/projects/JIA_PROcap/JIA_PROcap_mapping/seq_merged/denr_greater_than_1rpb_tx.bed.gz", 'rt') as f:
        for line in f:
            temp = line.strip().split("\t")
            key = temp[0]
            geneID = temp[3]
            TSS = int(temp[1])
            TES = int(temp[2])
            strand = temp[4]
            hashMapTSS[key].append([TSS, TES, 0, 0, strand, geneID])
    print(hashMapTSS.keys())

    #go through BigWig file and create hashmap where key is chromosome and value is a list
    #list has 6 values: [TSS, TES, # Pol II in TSS, # Pol II in TES, strand, geneID]

    print(f"Total number of lines: {line_count}")
    print("Header Information:")
    print(f"Chromosomes: {bw.chroms()}")
    print(f"Header Details: {bw.header()}")

    for chromosome in bw.chroms():
        if chromosome in hashMapTSS:
            for elem in hashMapTSS[chromosome]:
                try:
                    TSS = elem[0]
                    TES = elem[1]
                    strand = elem[4]
                    if strand == '+':
                        intervals = bw.intervals(chromosome, TSS - 150, TES - 300)
                    else:
                        continue
                    if intervals:
                        for interval in intervals:
                            if strand == '+':
                                if interval[0] <= elem[0] + 150:
                                    elem[2] += 1
                                else:
                                    elem[3] += 1
                except Exception as e:
                    print(f"Error fetching intervals: {e}")

    for chromosome in bw_neg.chroms():
        if chromosome in hashMapTSS:
            for elem in hashMapTSS[chromosome]:
                try:
                    TSS = elem[0]
                    TES = elem[1]
                    strand = elem[4]
                    if strand == '-':
                        intervals = bw_neg.intervals(chromosome, TSS + 300, TES + 150)
                    else:
                        continue
                    if intervals:
                        for interval in intervals:
                            if strand == '-':
                                if interval[0] >= TES - 150:
                                    elem[2] += 1
                                else:
                                    elem[3] += 1
                except Exception as e:
                    print(f"Error fetching intervals: {e}")
    output_file = f"/home2/npp8/data/seq{i}.pkl"
    with open(output_file, 'wb') as file:
        pickle.dump(hashMapTSS, file)

    print(f"hashMapTSS has been successfully saved to {output_file}")
    bw.close()


    input_file = f"/home2/npp8/data/seq{i}.pkl"
    output_csv = f"/home2/npp8/data/seq{i}_pause_index.csv"

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
                    denominator = (TES - 300) - (TSS + 150)
                    if denominator != 0:
                        ratio = (pol_ii_tss / 300) / (pol_ii_tes / denominator)
                        if ratio < 0:
                            print(pol_ii_tss)
                            print(pol_ii_tes)
                            print(denominator)
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

