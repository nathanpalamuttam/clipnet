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
#list has 4 values: [TSS, TES, # Pol II in TSS, # Pol II in TES]

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
                    intervals = bw.intervals(chromosome, TSS + 300, TES + 150)
                if intervals:
                    for interval in intervals:
                        if strand == '+':
                            if interval[0] <= elem[0] + 150:
                                elem[2] += 1
                            else:
                                elem[3] += 1
                        else:
                            if interval[0] >= TES - 150:
                                elem[2] += 1
                            else:
                                elem[3] += 1
            except Exception as e:
                print(f"Error fetching intervals: {e}")


output_file = "/home2/npp8/data/hashMapTSS.pkl"
with open(output_file, 'wb') as file:
    pickle.dump(hashMapTSS, file)

print(f"hashMapTSS has been successfully saved to {output_file}")
bw.close()
