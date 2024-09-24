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
    chrom_sizes = bw.chroms()
    for chromosome in bw.chroms():
        if chromosome in hashMapTSS:
            for elem in hashMapTSS[chromosome]:
                try:
                    promoter_start = max(TSS - 150, 0)  # Ensure start is not less than 0
                    promoter_end = min(TSS + 150, chrom_sizes[chromosome] )  # Ensure end is within chromosome length
                    gene_body_start = min(TSS + 300, chrom_sizes[chromosome] )  # Start of gene body
                    gene_body_end = max(TES - 300, 0)  # End of gene body, ensuring it's not negative

                    if strand == '+':
                        # Check if intervals are valid before querying stats
                        if promoter_start < promoter_end:
                            elem[2] = bw.stats(chromosome, promoter_start, promoter_end)[0]  # Get Pol II in TSS
                        if gene_body_start < gene_body_end:
                            elem[3] = bw.stats(chromosome, gene_body_start, gene_body_end)[0]
                    else:
                        continue
                # if intervals:
                #     for interval in intervals:
                #         if strand == '+':
                #             if interval[0] <= elem[0] + 150:
                #                 elem[2] += 1
                #             else:
                #                 elem[3] += 1
                except Exception as e:
                    print(f"Error fetching intervals: {e}")
    print("DONE")
    chrom_sizes_neg = bw_neg.chroms() 
    for chromosome in bw_neg.chroms():
        if chromosome in hashMapTSS:
            for elem in hashMapTSS[chromosome]:
                try:
                    TSS = elem[0]
                    TES = elem[1]
                    strand = elem[4]

                    # Define intervals with boundary checks for negative strand
                    promoter_start_neg = max(TES - 150, 0)  # Ensure start is not less than 0
                    promoter_end_neg = min(TES + 150, chrom_sizes_neg[chromosome] )  # Ensure end is within chromosome length
                    gene_body_start_neg = min(TSS + 300, chrom_sizes_neg[chromosome] )  # Start of gene body
                    gene_body_end_neg = max(TES - 300, 0)  # End of gene body, ensuring it's not negative

                    if strand == '-':
                        # Check if intervals are valid before querying stats
                        if promoter_start_neg < promoter_end_neg:
                            elem[2] = bw_neg.stats(chromosome, promoter_start_neg, promoter_end_neg)[0]  # Pol II in promoter
                        if gene_body_start_neg < gene_body_end_neg:
                            elem[3] = bw_neg.stats(chromosome, gene_body_start_neg, gene_body_end_neg)[0]  # Pol II in gene body

                    # if intervals:
                    #     for interval in intervals:
                    #         if strand == '-':
                    #             if interval[0] >= TES - 150:
                    #                 elem[2] += 1
                    #             else:
                    #                 elem[3] += 1
                except Exception as e:
                    print(f"Error fetching intervals: {e}")
    output_file = f"/home2/npp8/data/seq{i}_run2.pkl"
    with open(output_file, 'wb') as file:
        pickle.dump(hashMapTSS, file)

    print(f"hashMapTSS has been successfully saved to {output_file}")
    bw.close()


    input_file = f"/home2/npp8/data/seq{i}_run2.pkl"
    output_csv = f"/home2/npp8/data/seq{i}_pause_index_run2.csv"

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
                if pol_ii_tes != 0 and pol_ii_tss is not None and pol_ii_tes is not None:
                    
                    ratio = (pol_ii_tss) / (pol_ii_tes)
                else:
                    ratio = 'undefined'  # or some other handling for this specific case
                
                # if ratio == 'undefined' or ratio == 0:
                #     print(strand)
                
                # Write the row to the CSV file
                csv_writer.writerow([key, gene_id, ratio, strand, TSS, TES])

        print(f"Data successfully written to {output_csv}")

