### RASMUS HENRIKSEN JUNE 2020 ###
# Im using 3.6 BC09 as test file since its smaller than the rest..

# I chose to remove the contamination sources, to get a better overview.

samtools view -b -h BC09.aln_hg19.bam -L BC09_1000_cov.bed > BC09_hg19.bam

# I began to convert my NanoCircle2020.py script into argparse using the structure of the created arguments from april 2020

# However during the changes and additions in Simple_Cmd i found that some reads were still chimeric.

