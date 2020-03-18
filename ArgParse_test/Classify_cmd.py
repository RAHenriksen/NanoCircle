import argparse
import edlib as ed
import pandas as pd
import seaborn as sns
import numpy as np
import pybedtools as bt
from Bio import SeqIO
import re
import os
import pysam as ps
import matplotlib.pyplot as plt
import pylab as plt2
import time
from collections import defaultdict
from matplotlib.patches import Rectangle
from collections import Counter
import pybedtools
import math


class utils:
    def __init__(self,cigar_str):
        self.cigar_str = cigar_str

    def CIGAR_len(self):
        """
        :param cigar_str: Cigar string from alignment information from SAM/BAM format
        :return: Length of the alignment from CIGAR string
        """
        count = 0
        for i in re.findall(r'(\d+)([A-Z]{1})', str(self.cigar_str)):
            #Match, Insertion, Hardclip
            if i[-1] == "M" or i[-1] == "I" or i[-1] == "H":
                count += int(i[0])
            #Deletion
            elif i[-1] == "D":
                count -= int(i[0])
        return count


class Read_Filter:
    """ Class which extract alignment position information for filtered and unfiltered reads"""
    def __init__(self, bamfile,tag,MapQ):
        self.bamfile = bamfile
        self.tag = tag
        self.MapQ = MapQ

    def SamBamParser(self):
        """ Loads the sam / bam file """
        #add .fetch(reg,start,end,multiple_iterators=True)
        file = ps.AlignmentFile(self.bamfile, "rb")
        return file

    def Filter(self):
        """
        Filters the reads, based on the cigar operation, read tag and read quality given to the class
        :return: List with Pysam objects of the reads
        """
        file=self.SamBamParser()
        Simple_reads = []
        Chimeric_reads = []
        for read in file.fetch():
            print("read",read)
            # filtering on the cigar string operation, S = 4, H = 5
            if read.cigartuples[0][0] == 4 or read.cigartuples[-1][0] == 4:
                # filtering on mapping quality
                if read.mapping_quality >= self.MapQ:
                    # filtering on those reads with a pointer to part of read with tag.
                    # With prim pointing to supp, and supp pointing back to prim with "SA" tag
                    if read.has_tag(self.tag):
                        Tag_full = read.get_tag(self.tag).split(';')[:-1]
                        Tag_elem = [elem.split(',') for elem in Tag_full]

                        if len(Tag_elem) == 1 and read.reference_name == Tag_elem[0][0]:
                            Simple_reads.append(read)
                        else:
                            # if the SA tag have several alignments, it takes the unique set #set()
                            # of every first element of the tag elem  list(zip(*Tag_elem))[0]
                            chr = set(list(zip(*Tag_elem))[0])

                            # if the len of the unique is 1 we assume its a simple with multiple SA
                                # think about chimeric which might just be two regions from same chromomsome
                                    # for another time
                            if len(chr) == 1:
                                Simple_reads.append(read)

                            #Otherwise chimeric
                            else:
                                Chimeric_reads.append(read)

        return Simple_reads,Chimeric_reads
    #    if self.Method == "Simple":
    #        return Simple_reads
    #    elif self.Method == "Chimeric":
    #        return Chimeric_reads

if __name__ == '__main__':
    #print(Merged_dict)

    Read_class=Read_Filter("/isdata/common/wql443/NanoCircle/BC02.aln_hg19.bam","SA",30)
    simple,chimeric = Read_class.Filter()
    print([i.query_name for i in simple])
    """print(simple[0].reference_name)
    print(simple[0].query_name)
    print(simple[0].get_tag("SA"))
    print("works")
    print(chimeric[0].reference_name)
    print(chimeric[0].query_name)
    print(chimeric[0].get_tag("SA"))"""