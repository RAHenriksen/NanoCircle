#!/usr/bin/env python
import edlib as ed
import pandas as pd
import seaborn as sns
import numpy as np
from Bio import SeqIO
import re
import os
import pysam as ps
import matplotlib.pyplot as plt
import pylab as plt2
import time
from collections import defaultdict
from matplotlib.patches import Rectangle
import pybedtools
#from pybedtools import BedTool

def density_plot(dict1):
    """ plots the distribution of split-reads across a chromosome region """
    x1 = []
    x2 = []
    #append the positions to list
    for k,v in dict1.items():
        x1.append(v[0])
        x2.append(v[1])
    colors=["dodgerblue","darkorange"]
    labels = ["Primary alignment","Supplementary alignment"]
    plt.hist([x1,x2],bins=20,density=True,color=colors,label=labels)
    #plt.title(str(len(Read_list)) + " Soft-clipped reads " + label)
    plt.legend(loc='upper center')
    plt.show()

def Softclipped_NoSA_count(bamfile,tag_str):
    """Returns list with length of soft_clipped reads with and without SA"""
    SA_list = []
    No_list = []
    for read in bamfile.fetch():
        # the left most softclip read
        if read.cigar[0][0] == 4:  #cigar 4 equals Softclipped in Sam format
            if read.has_tag(tag_str):
                SA_list.append(read.cigar[0][1])
            if not read.has_tag(tag_str):
                No_list.append(read.cigar[0][1])
        # the right most split read
        if read.cigar[-1][0] == 4:
            if read.has_tag(tag_str):
                SA_list.append(read.cigar[-1][1])
            if not read.has_tag(tag_str):
                No_list.append(read.cigar[-1][1])
    return SA_list, No_list

def Read_dist(Read_list,bins,label,figlabel,int1,int2):
    plt2.hist(Read_list, bins=bins, alpha=0.5, label=label, color="blue")
    plt2.hist([], color='white',
              label= "{:.2f}".format(sum(i > int1 for i in Read_list) / len(Read_list) * 100)
                     + "%" + " of soft-clipped part w. length above %s bp" % str(int1))
    plt2.hist([], color = 'white',
              label="{:.2f}".format(sum(i < int2 for i in Read_list) / len(Read_list) * 100)
                    + "%" + " of soft-clipped part w. length below %s bp" % str(int2))
    plt2.title(str(len(Read_list)) + " Soft-clipped reads " + label)
    plt2.legend()
    plt2.xlabel("Soft-clipped length")
    plt2.ylabel("Read count")
    plt2.savefig(figlabel+".jpg")
    plt2.show()

### PLOTTING CODE ##
#test = Coord_dict(bamfile, False, "chr1", 243928620, 243938331)
#density_plot(test)
#test = Coord_dict(bamfile, True, "chr1", 243928620, 243938331)
#density_plot(test)

#bamfile2 = ps.AlignmentFile("chr1_243928620_243938331_region.bam", "rb")

#SA_list,NoSA_list = Softclipped_NoSA_count(bamfile2,"SA")

#Read_dist(SA_list,100,"with supplementary alignment","SA",1000,100)
#Read_dist(NoSA_list,100,"without supplementary alignment","NoSA",1000,100)

print("dette lort virker efter mounting shit")
