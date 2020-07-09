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
    def __init__(self, bamfile,dir,MapQ):
        self.bamfile = bamfile
        self.dir = dir
        self.MapQ = MapQ

    def Filter(self):
        """
        Filters the reads, based on the cigar operation, read tag and read quality given to the class
        :return: List with Pysam objects of the reads
        """
        file = Utils.SamBamParser(self.bamfile)

        #creates alignment files, with the different reads seperated
        ps_simple = ps.AlignmentFile("{0}/Simple_reads.bam".format(str(self.dir)), "wb", template=file)
        ps_chimeric = ps.AlignmentFile("{0}/Chimeric_reads.bam".format(str(self.dir)), "wb", template=file)

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

                        Prim_pos = ps.AlignedSegment.get_reference_positions(read)

                        #print(" TAG ELEM",read.query_name,read.reference_name,Prim_pos[0],Prim_pos[-1],Tag_elem)

                        # one supp alignment
                        if len(Tag_elem) == 1:
                            # ensures the supp is aligning to same chromosome
                            if read.reference_name == str(Tag_elem[0][0]):
                                # Handle the reads with dist above 100 k as chim with single supp, but it could in theory be simple
                                if abs(Prim_pos[0]-int(Tag_elem[0][1])) > 100000 and \
                                        abs(Prim_pos[1]-(int(Tag_elem[0][1])+Utils.CIGAR_len(Tag_elem[0][3]))) > 100000:

                                    #ensures the chimeric read supp alignment is also above mapQ
                                    if int(Tag_elem[0][4]) >= self.MapQ:
                                        ps_chimeric.write(read)
                                    #if not we discard the read
                                    else:
                                        continue

                                else:
                                    ps_simple.write(read)

                            #single supp to another chr
                            else:
                                ps_chimeric.write(read)

                        # multiple supp align
                        elif len(Tag_elem) > 1:
                            #print(" Chim tag", read.query_name,read.reference_name, Tag_elem)

                            # if the SA tag have several alignments, it takes the unique set #set() of
                            # the given elements with zip
                            Prim_interval = [Prim_pos[0] - 50000, Prim_pos[1] + 50000]

                            chr = list(set(list(zip(*Tag_elem))[0]))
                            pos_start = list(list(zip(*Tag_elem))[1])
                            #maps cigar list on all cigar strings found with zip*
                            read_len = list(map(Utils.CIGAR_len,list(list(zip(*Tag_elem))[3])))
                            pos_interval = []

                            for i in range(len(pos_start)):
                                pos_interval.append(int(pos_start[i]))
                                pos_interval.append(int(pos_start[i]) + read_len[i])

                            # if the len of the unique set of chr is 1 its the same chromosome for all the SA
                            if len(chr) == 1:
                                # if all supp is the same chromosome as prim
                                if read.reference_name == str(chr[0]):
                                    #due to the distance, if some supp alignment are outside the +- 50 k interval of prim its chimeric
                                    if any(e < Prim_interval[0] or Prim_interval[1] < e for e in pos_interval):
                                        ps_chimeric.write(read)
                                    else:
                                        ps_simple.write(read)

                                # if all supp is to another chrom
                                else:
                                    ps_chimeric.write(read)

                            #multiple different chr
                            else:
                                ps_chimeric.write(read)

        ps_simple.close()
        ps_chimeric.close()




if __name__ == '__main__':
    #print(Merged_dict)

    Read_class=Read_Filter("/isdata/common/wql443/NanoCircle/ArgParse_test/BC10.aln_hg19.bam","temp_reads",30)
    Read_class.Filter()
    """

    prim_pos = [45788715,45790577]
    prim_interval = [prim_pos[0]-50000,prim_pos[1]+50000]
    Supp_start = [45788851, 45788705, 45788706, 54265141, 54265145, 54265621, 54265141]
    Supp_pos = []

    print(any(prim_interval[0] < e < prim_interval[1] for e in Supp_start))
    print(any(e < prim_interval[0] or prim_interval[1] < e for e in Supp_start))

    a = ['45788851', '45788705', '45788706', '54265141', '54265145', '54265621', '54265141']
    b = [1676, 1762, 1134, 960, 974, 524, 323]

    c= []
    for i in range(len(a)):
        c.append(int(a[i]))
        c.append(int(a[i])+b[i])
    print(c)
    a = [ i for i in a[i]]"""