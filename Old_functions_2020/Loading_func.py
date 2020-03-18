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


def Genome_Cov_identification(bamfile,bam_name,overlap,mapQ):
    Bedfile = pybedtools.example_bedtool(str(os.getcwd()) +"/"+bam_name)
    Cov = Bedfile.genome_coverage(bg=True)
    Merged = Cov.merge(d=overlap)

    for region in Merged:
        Chr = region[0]
        start = region[1]
        end = region[2]
        if len(Single_coord(bamfile, mapQ,str(Chr), int(start), int(end))) == 0:
            continue
        else:
            print(Chr + ":" + start + "-" + end)
            print(Single_coord(bamfile, mapQ, str(Chr), int(start), int(end)))

def File_identification(file_name,mapQ):
    count = 0
    with open(file_name) as f:
        for line in f:
            line_value = line.strip().split()
            coord = line_value[0]
            start = line_value[1]
            end = line_value[2]
            cov = line_value[3]
            coord_dict = Single_coord(bamfile,mapQ, str(coord), int(start), int(end))

            if bool(coord_dict) == True:
                print(coord+":"+start+"-"+end)
                print("dictionary",coord_dict)
                count += 1
            else:
                #those dict empy, possible due to complex circle formation
                continue
    print("the number of circles is",count)

Genome_Cov_identification(bamfile,"BC05.aln_hg19.bam",500,60)

File_identification("BC05.ge_mean5.bdg",60)
