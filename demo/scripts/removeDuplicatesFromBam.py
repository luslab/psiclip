# Remove PCR duplicates from BAM file
# Based on rbc and start position
# Charlotte Capitanchik
# Concept:
# Parse sam file, make a tuple representing each unique cDNA - (rbc, start, strand)
# this becomes a key in a dictionary where the values are the read names
# Then to make the de-deduplicated sam file take the first item in each 
# value list and retrieve those reads.

import pysam
import sys
from collections import defaultdict

bamfile = sys.argv[1] 
#"/Users/capitac/Desktop/m1a_HEK293_beta_1.Aligned.out.sorted.bam"
outfile = sys.argv[2] 
#"/Users/capitac/Desktop/m1a_HEK293_beta_1.Aligned.out.sorted.deduplicated.bam"

samfile = pysam.AlignmentFile(bamfile,"rb")
read_dict = defaultdict(list)

for read in samfile.fetch():
    if read.is_unmapped == False:
        rbc = read.query_name.split(":")[-1]
        start = "".join(list(map(str,read.get_reference_positions())))
        strand = read.is_reverse
        read_dict[(rbc,start,strand)].append(read.query_name)
    else:
        continue

unique_reads = [value[0] for value in read_dict.values()]
unique_reads = set(unique_reads)

dedup_bam = pysam.AlignmentFile(outfile, "wb", template=samfile)

for read in samfile.fetch():
    if read.query_name in unique_reads:
        dedup_bam.write(read)
