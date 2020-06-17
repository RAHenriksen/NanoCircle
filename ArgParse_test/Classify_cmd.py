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

import Utils

class Read_Filter:
    """ Class which extract alignment position information for filtered and unfiltered reads"""
    def __init__(self, bamfile,MapQ):
        self.bamfile = bamfile
        self.MapQ = MapQ

    def Filter(self):
        """
        Filters the reads, based on the cigar operation, read tag and read quality given to the class
        :return: List with Pysam objects of the reads
        """
        file = Utils.SamBamParser(self.bamfile)

        #creates alignment files, with the different reads seperated
        ps_simple = ps.AlignmentFile("temp_reads/Simple_reads.bam", "wb", template=file)
        ps_chimeric = ps.AlignmentFile("temp_reads/Chimeric_reads.bam", "wb", template=file)

        for read in file.fetch():
            # filtering on the cigar string operation, S = 4
            if read.cigartuples[0][0] == 4 or read.cigartuples[-1][0] == 4:
                # filtering on mapping quality
                if read.mapping_quality >= self.MapQ:
                    # filtering on those reads with a pointer to part of read with tag.
                    # With prim pointing to supp, and supp pointing back to prim with "SA" tag
                    if read.has_tag("SA"):
                        Tag_full = read.get_tag("SA").split(';')[:-1]
                        # creating a nested list, with each tag elem as list elem
                        Tag_elem = [elem.split(',') for elem in Tag_full]

                        if len(Tag_elem) == 1:
                            # ensures the supp is aligning to same chromosome
                            if read.reference_name == str(Tag_elem[0][0]):
                                print("SIMPLE", read.reference_name, Tag_elem[0][0],
                                      Tag_elem)
                                ps_simple.write(read)
                            elif read.reference_name != str(Tag_elem[0][0]):
                                ps_chimeric.write(read)

                        #multiple supp align
                        elif len(Tag_elem) > 1:
                            # if the SA tag have several alignments, it takes the unique set #set()
                            # of every first element of the tag elem  list(zip(*Tag_elem))[0]
                            chr = list(set(list(zip(*Tag_elem))[0]))
                            # if the len of the unique is 1 we assume its a simple with multiple SA
                                # think about chimeric which might just be two regions from same chromomsome
                                    # for another time
                            #print(chr)
                            #print(len(chr))
                            if len(chr) == 1 and str(chr[0]) == read.reference_name:
                                ps_simple.write(read)
                                #print("SIMPLE", read.reference_name, Tag_elem)
                            #Otherwise chimeric
                            else:
                                ps_chimeric.write(read)
                                #Chimeric_reads.append(read)

        ps_simple.close()
        ps_chimeric.close()


if __name__ == '__main__':
    #print(Merged_dict)

    Read_class=Read_Filter("/isdata/common/wql443/NanoCircle/BC02.aln_hg19.bam",30)
    Read_class.Filter()
